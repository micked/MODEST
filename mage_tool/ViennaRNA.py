#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Python wrapper for the Vienna RNA Package by Andreas R. Gruber, Ronny Lorenz,
Stephan H. Bernhart, Richard Neubock, and Ivo L. Hofacker (NAR, 2008).
"""

import re
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

    RNAFOLD    = "RNAfold"
    RNACOFOLD  = "RNAcofold"
    RNASUBOPT  = "RNAsubopt"
    RNAEVAL    = "RNAeval"

    def __init__(self, version="auto"):
        """TODO"""

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
            if len(self.versions) > 1:
                self.fallbackversion = self.versions[1]
            else:
                self.fallbackversion = False
        else:
            if not version in self.versions:
                raise Exception("Specified ViennaRNA version {} is not availble. Found: {}".format(version, self.versions))
            self.version = version

    def fold(self, seq, noPS=True, **kwargs):
        """RNAfold"""
        version = self.version

        if self.version == "native":
            if not noPS or kwargs:
                version = self.fallbackversion
                if not self.fallbackversion:
                    raise Exception("ViennaRNA native cannot take any arguments, and no fallback available")
            else:
                return self.RNA.fold(seq)

        kwargs["noPS"] = noPS

        options = self.versionize_kwargs(version, **kwargs)

        output = self.run_command([self.RNAFOLD] + options, seq)
        output = output.splitlines()
        if len(output) > 1:
            return self.output_to_brackets_and_dG(output[1])
        else:
            print seq
            print " ".join([self.RNAFOLD] + options)
            return None, None

    def cofold(self, seq, noPS=True, **kwargs):
        """RNAcofold"""
        version = self.version

        if self.version == "native":
            if not noPS or kwargs:
                version = self.fallbackversion
                if not self.fallbackversion:
                    raise Exception("ViennaRNA native cannot take any arguments, and no fallback available")
            else:
                return self.RNA.cofold(seq)

        kwargs["noPS"] = noPS

        options = self.versionize_kwargs(version, **kwargs)

        output = self.run_command([self.RNACOFOLD] + options, seq)
        output = output.splitlines()

        if len(output) > 1:
            return self.output_to_brackets_and_dG(output[1])
        else:
            return None, None

    def subopt(self, seq1, seq2, **kwargs):
        """RNAsubopt"""
        version = self.version

        if self.version == "native":
            raise NotImplementedError("No native subopt yet")

        options = self.versionize_kwargs(version, **kwargs)

        seq = "&".join([seq1, seq2])
        output = self.run_command([self.RNASUBOPT] + options, seq)
        output = output.splitlines()

        if len(output) > 1:
            return [self.output_to_brackets_and_dG(line) for line in output[1:]]
        else:
            return []

    
    def energy(self, seq, brackets, **kwargs):
        """RNAsubopt"""
        version = self.version

        if self.version == "native":
            raise NotImplementedError("No native energy yet")

        options = self.versionize_kwargs(version, **kwargs)

        inp = "\n".join([seq, brackets])
        output = self.run_command([self.RNAEVAL] + options, inp)
        output = output.splitlines()

        if len(output) > 1:
            return self.output_to_brackets_and_dG(output[1])[1]
        else:
            return None


    def versionize_kwargs(self, version, **kwargs):
        """Turn kwargs into a list of options fitting current version"""
        options = []
        if version == 2:
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


if __name__ == "__main__":
    import doctest
    doctest.testmod()