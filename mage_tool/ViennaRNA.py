#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Python wrapper for the Vienna RNA Package by Andreas R. Gruber, Ronny Lorenz,
Stephan H. Bernhart, Richard Neubock, and Ivo L. Hofacker (NAR, 2008).
"""

from __future__ import print_function

import re
import os
from subprocess import Popen, PIPE


def mfe(seq):
    """Wrapper around RNAcofold, return only dG"""
    return ViennaRNA().cofold(seq)[1]


def brackets_to_basepairing(brackets):
    """Convert a bracket string notation to numbered pairs

    Returned is a list of length of strands, a list of x values and a list of
    y_values:

        >>> brackets_to_basepairing("((...))((.(....&.)))")
        ([15, 4], [0, 1, 7, 8, 10], [6, 5, 18, 17, 16])

    The function will test for missing brackets:

        >>> brackets_to_basepairing("((...))((.(....&.))")
        Traceback (most recent call last):
         ...
        ValueError: Missing "y" bracket in ((...))((.(....&.))

    """
    bp_x = list()
    last = list()
    bp_y = [0] * brackets.count(")")
    strands = [len(s) for s in brackets.split("&")]

    for i,letter in enumerate(brackets.replace("&","")):
        if letter == "(":
            bp_x.append(i)
            last.append(i)
        elif letter == ")":
            try:
                matching_bracket = last.pop()
                y_loc = bp_x.index(matching_bracket)
                bp_y[y_loc] = i
            except IndexError:
                a = "x"
                if brackets.count("(") > brackets.count(")"): a ="y"
                raise ValueError("Missing \"{}\" bracket in {}".format(a, brackets))

    return strands, bp_x, bp_y


def basepairing_to_brackets(bp_x, bp_y, seq=None, strands=None):
    """Convert two lists of matching basepair lists to bracket notation.

    Either a strands-list must be supplied:

        >>> strands, x, y = ([15, 4], [0, 1, 7, 8, 10], [6, 5, 18, 17, 16])
        >>> basepairing_to_brackets(x, y, strands=strands)
        '((...))((.(....&.)))'

    Or the original sequence:

        >>> strands, x, y = brackets_to_basepairing(".(((....)))..")
        >>> seq = "AUCGACUACGACU"
        >>> basepairing_to_brackets(x, y, seq=seq)
        '.(((....)))..'

    """
    if strands is not None:
        cut_points = [strands[i] + (strands[i-1] if i else 0) for i in range(len(strands)-1)]
        l = sum(strands)
    elif seq is not None:
        cut_points = [m.start() for m in re.finditer(r"\&", seq)]
        l = len(seq) - len(cut_points)
    else:
        raise AttributeError("Either seq or strands must be supplied.")

    brackets = ["."] * l
    for x, y in zip(bp_x, bp_y):
        brackets[x] = "("
        brackets[y] = ")"

    [brackets.insert(m, "&") for m in cut_points]

    return "".join(brackets)


class ViennaRNA:

    def __init__(self, version="auto", prefix=""):
        """TODO"""

        self.RNAFOLD    = os.path.join(prefix, "RNAfold")
        self.RNACOFOLD  = os.path.join(prefix, "RNAcofold")
        self.RNASUBOPT  = os.path.join(prefix, "RNAsubopt")
        self.RNAEVAL    = os.path.join(prefix, "RNAeval")
        self.RNADUPLEX  = os.path.join(prefix, "RNAduplex")

        self.versions = list()
        #Try native version
        try:
            import RNA
            structure, dG = RNA.fold("CCCCCAAAGGGGG")
            #Extra output check
            if len(structure) != 13:
                raise Exception("Invalid RNAfold output")

            self.RNA = RNA
            self.versions.append("native")
        except Exception:
            pass

        #Try version 2.X
        try:
            self.run_command([self.RNAFOLD, "--version"])
            self.versions.append(2)
        except Exception:
            pass

        #Try version 1.X
        try:
            self.run_command([self.RNAFOLD, "-noPS"], "CCCCCAAAGGGGG")
            self.versions.append(1)
        except Exception:
            pass

        if not self.versions:
            raise Exception("No valid ViennaRNA installation found.")

        if version == "auto":
            self.version = self.versions[0]
        elif version == "nonnative":
            if 2 in self.versions:
                self.version = 2
            elif 1 in self.versions:
                self.version = 1
            else:
                raise ValueError("No nonnative ViennaRNA found.")
        else:
            if not version in self.versions:
                raise ValueError("Specified ViennaRNA version {} is not availble. Found: {}".format(version, self.versions))
            self.version = version

    def fold(self, seq, d=1):
        """RNAfold"""

        if self.version == "native":
            self.RNA.cvar.dangles = d
            self.RNA.cvar.cut_point = -1
            f = self.RNA.fold(seq)
            self.RNA.free_arrays()
            return f

        kwargs = {"noPS": True, "d": d}
        options = self.versionize_kwargs(kwargs)

        output = self.run_command([self.RNAFOLD] + options, seq)
        output = output.splitlines()
        if len(output) > 1:
            return self.output_to_brackets_and_dG(output[1])
        else:
            # print seq
            # print " ".join([self.RNAFOLD] + options)
            return None, None

    def cofold(self, seqs, d=1):
        """RNAcofold"""
        if type(seqs) is list:
            seqs = "&".join(seqs)

        cps = seqs.count("&")
        if cps > 1:
            raise ValueError("Use maximum 1 cut-point in cofold.")
        elif cps == 0:
            return self.fold(seqs, d=d)

        if self.version == "native":
            self.RNA.cvar.dangles   = d
            #Find cut point and set
            cp = seqs.find("&")
            self.RNA.cvar.cut_point = cp + 1
            #Erase cut point marker
            seqs = seqs.replace("&", "")
            #Fold
            f = self.RNA.cofold(seqs)
            #Add cut point back
            f[0] = self.add_cut_point(f[0], cp)
            self.RNA.cvar.cut_point = -1
            return f

        kwargs = {"noPS": True, "d": d}

        options = self.versionize_kwargs(kwargs)

        output = self.run_command([self.RNACOFOLD] + options, seqs)
        output = output.splitlines()

        if len(output) > 1:
            return self.output_to_brackets_and_dG(output[1])
        else:
            return None, None


    def duplex(self, seqs, d=1):
        """RNAcofold"""
        try:
            seqs = seqs.split("&")
        except AttributeError:
            pass

        if len(seqs) != 2:
            raise Exception('Duplexfold can only fold 2 sequences')

        if self.version == "native":
            self.RNA.cvar.dangles   = d
            f = self.RNA.duplexfold(*seqs)
            return f.structure,f.energy

        kwargs = {"noPS": True, "d": d}

        options = self.versionize_kwargs(kwargs)

        output = self.run_command([self.RNADUPLEX] + options, '\n'.join(seqs))
        output = output.splitlines()

        if len(output) > 1:
            return self.output_to_brackets_and_dG(output[1])
        else:
            return None, None

    def subopt(self, seqs, d=1, e=1):
        """RNAsubopt"""
        if type(seqs) is list:
            seqs = "&".join(seqs)

        cps = seqs.count("&")
        if cps > 1:
            raise ValueError("Use maximum 1 cut-point in subopt.")

        if self.version == "native":
            mfe = self.cofold(seqs, d=d)[1]

            cp = seqs.find("&")
            self.RNA.cvar.cut_point = -1
            if cps:
                self.RNA.cvar.cut_point = cp + 1
                seqs = seqs.replace("&", "")

            self.RNA.cvar.dangles = d
            sub = self.RNA.subopt(seqs, None, int(round(e*100)))

            output = list()
            for i in range(sub.size()):
                if cps:
                    brackets = self.add_cut_point(sub.get(i).structure, cp)
                else:
                    brackets = sub.get(i).structure

                energy = sub.get(i).energy

                #print mfe, e, energy

                if mfe+e > energy:
                    output.append((brackets, energy))

            return output

        kwargs = {"d": d, "e": e}
        options = self.versionize_kwargs(kwargs)

        output = self.run_command([self.RNASUBOPT] + options, seqs)
        output = output.splitlines()

        if len(output) > 1:
            return [self.output_to_brackets_and_dG(line) for line in output[1:]]
        else:
            return []

    def energy(self, seq, brackets, d=1):
        """RNAeval"""
        if type(seq) is list:
            seq = "&".join(seq)

        if len(seq) != len(brackets):
            raise ValueError("Sequence and structure have unequal length.")

        if self.version == "native":
            cp = seq.find("&")
            self.RNA.cvar.cut_point = -1
            if cp != -1:
                self.RNA.cvar.cut_point = cp + 1

            self.RNA.cvar.dangles = d
            output = self.RNA.energy_of_struct(seq, brackets)
            self.RNA.cvar.cut_point = -1
            return output

        kwargs = {"d": d}
        options = self.versionize_kwargs(kwargs)

        inp = "\n".join([seq, brackets])
        output = self.run_command([self.RNAEVAL] + options, inp)
        output = output.splitlines()

        if len(output) > 1:
            return self.output_to_brackets_and_dG(output[1])[1]
        else:
            return None

    def add_cut_point(self, structure, cp):
        """Add a cut point (&) to the structure"""
        return structure[:cp] + "&" + structure[cp:]

    def versionize_kwargs(self, kwargs):
        """Turn kwargs into a list of options fitting current version"""
        options = []
        if self.version == 2:
            for kw in kwargs:
                if len(kw) == 1:
                    options.append("-" + kw)
                    if not (kwargs[kw] is True or kwargs[kw] is False):
                        options[-1] += str(kwargs[kw])
                else:
                    options.append("--" + kw)
                    if not (kwargs[kw] is True or kwargs[kw] is False):
                        options.append(str(kwargs[kw]))
        else:
            for kw in kwargs:
                options.append("-" + kw)
                if not (kwargs[kw] is True or kwargs[kw] is False):
                    if kw in ["d"]:
                        options[-1] += str(kwargs[kw])
                    else:
                        options.append(str(kwargs[kw]))

        return options

    def run_command(self, cmd, input_string=""):
        """Run the specified command and return output"""
        p = Popen(cmd, stdout=PIPE, stdin=PIPE, stderr=PIPE, universal_newlines=True)
        out, stderr = p.communicate(input=input_string)
        if p.returncode:
            raise Exception("Cmd {} failed: {}".format(cmd[0], stderr))
        return out

    def output_to_brackets_and_dG(self, outputline):
        """Single output line to brackets and dG"""
        outputline = outputline.split()
        brackets = outputline[0]
        dG = float("".join(outputline[1:]).strip("()"))
        return brackets, dG


"""
Unittests
~~~~~~~~~
"""

import unittest


class InterfaceTests(unittest.TestCase):
    def setUp(self):
        self.rna = ViennaRNA()

    def test_fold(self):
        fold, mfe = self.rna.fold("ATCGATCGATCGATCG")
        self.assertIsInstance(mfe, float)
        self.assertEqual(fold.count("("), fold.count(")"))
        self.assertRegexpMatches(fold, "[\.\(\)]+")

    def test_cofold(self):
        fold, mfe = self.rna.cofold("GCATGCAT&ATGCATGC")
        self.assertEqual(fold, "((((((((&))))))))")


if __name__ == "__main__":
    import doctest
    print("Running doctests..")
    doctest.testmod()

    print()
    print("Running unittests..")
    v = ViennaRNA().version
    print("Testing version:", v)

    unittest.main()
