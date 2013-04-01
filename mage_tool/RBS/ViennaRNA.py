import re
from subprocess import Popen, PIPE

RNACOFOLD = "RNAcofold"

def run_command(cmd, input_string):
    """Run the specified command and return output"""
    p = Popen(cmd, stdout=PIPE, stdin=PIPE, stderr=PIPE, universal_newlines=True)
    out, stderr = p.communicate(input=input_string)
    if p.returncode:
        raise Exception("Cmd {} failed: {}".format(cmd[0], stderr))
    return out

try:
    run_command(["RNAcofold", "-noPS"], "")
except Exception:
    noPS = "--noPS"
else:
    noPS = "-noPS"

def rnacofold(seq, options=[noPS], calc_basepairing=True):
    """Run RNAcofold"""
    cmd = [RNACOFOLD] + options
    output = run_command(cmd, seq)
    brackets, dG = output_to_brackets_and_dG(output.splitlines()[1])
    
    out = {"dG": dG, "brackets": brackets}

    if calc_basepairing:
        out["strands"], out["bp_x"], out["bp_y"] = brackets_to_basepairing(brackets)

    return out


def mfe(seq):
    """Wrapper around RNAcofold, return only dG"""
    return rnacofold(seq, calc_basepairing=False)["dG"]


def output_to_brackets_and_dG(outputline):
    """Single output line to brackets and dG"""
    outputline = outputline.split()
    brackets = outputline[0]
    dG = float("".join(outputline[1:]).strip("()"))
    return brackets, dG


def brackets_to_basepairing(brackets):
    bp_x = list()
    last = list()
    bp_y = [0] * brackets.count(")")
    strands = [len(s) for s in brackets.split("&")]

    for i,letter in enumerate(brackets.replace("&",""), 1):
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
                raise Exception("Missing \"{}\" bracket in {}".format(a, brackets))

    return strands, bp_x, bp_y


"""
Modfied from Ribosome-Binding-Site-Calculator:
<https://github.com/hsalis/ribosome-binding-site-calculator>
<https://salis.psu.edu/software/>

Python wrapper for the Vienna RNA Package by Andreas R. Gruber, Ronny Lorenz, Stephan H. Bernhart, Richard Neubock,
and Ivo L. Hofacker (NAR, 2008).

This file is part of the Ribosome Binding Site Calculator.

The Ribosome Binding Site Calculator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The Ribosome Binding Site Calculator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Ribosome Binding Site Calculator.  If not, see <http://www.gnu.org/licenses/>.

This Python wrapper is written by Howard Salis. Copyright 2008-2009 is owned by the University of California Regents.
All rights reserved. :) Use at your own risk.
"""


debug = 0

#Class that encapsulates all of the functions from NuPACK 2.0
class ViennaRNA(dict):

    RT = 0.61597 #gas constant times 310 Kelvin (in units of kcal/mol)

    def __init__(self, sequence_list, material="rna37"):
        self.ran = 0
        nuc_regex = re.compile(r"^[ATGCU]+$", re.IGNORECASE)

        for seq in sequence_list:
            if not nuc_regex.match(seq):
                error_string = "Invalid letters found in inputted sequences. Only ATGCU allowed. \n Sequence is \"" + str(seq) + "\"."
                raise ValueError(error_string)

        self["sequences"] = sequence_list
        self["material"] = material


    def mfe(self, strands, Temp=37.0, dangles="all", outputPS=False):
        self["mfe_composition"] = strands

        if Temp <= 0: raise ValueError("The specified temperature must be greater than zero.")

        input_string = "&".join(self["sequences"])

        cmd = ["RNAcofold"]

        if not outputPS:
            cmd.append("--noPS")

        if dangles is "none":
            cmd.append("-d0")
        elif dangles is "some":
            cmd.append("-d1")
        elif dangles is "all":
            cmd.append("-d2")

        #Call ViennaRNA C programs
        output = self._run_command(cmd, input_string)

        if debug:
            print input_string
            print " ".join(cmd)
            print output

        output = output.split("\n")

        words = output[1].split()
        bracket_string = words[0]
        (strands,bp_x, bp_y) = self.convert_bracket_to_numbered_pairs(bracket_string)

        energy = float(words[len(words)-1].replace(")","").replace("(","").replace("\n",""))

        self["program"] = "mfe"
        self["mfe_basepairing_x"] = [bp_x]
        self["mfe_basepairing_y"] = [bp_y]
        self["mfe_energy"] = [energy]
        self["totalnt"]=strands


    def subopt(self, strands, energy_gap, Temp = 37.0, dangles="all", outputPS=False):
        self["subopt_composition"] = strands

        if Temp <= 0: raise ValueError("The specified temperature must be greater than zero.")

        input_string = "&".join(self["sequences"])

        cmd = ["RNAsubopt"]
        cmd.extend(["-e", str(energy_gap)])

        if dangles is "none":
            cmd.append("-d0")
        elif dangles is "some":
            cmd.append("-d1")
        elif dangles is "all":
            cmd.append("-d2")

        #Call ViennaRNA C programs
        output = self._run_command(cmd, input_string)

        if debug or outputPS:
            print input_string
            print " ".join(cmd)
            print output

        output = output.split("\n")

        self["subopt_basepairing_x"] = []
        self["subopt_basepairing_y"] = []
        self["subopt_energy"] = []
        self["totalnt"]=[]
        self["brackets"]=[]
        counter=0

        for line in output[1:]:
            if len(line) > 0:
                counter+=1
                words = line.split(" ")
                bracket_string = words[0]
                energy = float(words[len(words)-1].replace("\n",""))

                (strands,bp_x, bp_y) = self.convert_bracket_to_numbered_pairs(bracket_string)

                self["subopt_energy"].append(energy)
                self["subopt_basepairing_x"].append(bp_x)
                self["subopt_basepairing_y"].append(bp_y)
                self["brackets"].append(bracket_string)

        self["subopt_NumStructs"] = counter
        self["program"] = "subopt"


    def energy(self, strands, base_pairing_x, base_pairing_y, Temp = 37.0, dangles = "all"):
        self["energy_composition"] = strands

        if Temp <= 0: raise ValueError("The specified temperature must be greater than zero.")

        seq_string = "&".join(self["sequences"])
        strands = [len(seq) for seq in self["sequences"]]
        bracket_string = self.convert_numbered_pairs_to_bracket(strands,base_pairing_x,base_pairing_y)
        input_string = seq_string + "\n" + bracket_string + "\n"

        #print seq_string
        # print "te", bracket_string

        cmd = ["RNAeval"]

        #Set arguments
        if dangles is "none":
            cmd.append("-d0")
        elif dangles is "some":
            cmd.append("-d1")
        elif dangles is "all":
            cmd.append("-d2")

        #Call ViennaRNA C programs
        output = self._run_command(cmd, input_string)

        if debug:
            print input_string
            print " ".join(cmd)
            print output

        output = output.split("\n")

        self["energy_energy"] = []

        #Skip the unnecessary output lines
        words = output[1].split()
        energy = float(words[len(words)-1].replace("(","").replace(")","").replace("\n",""))

        self["program"] = "energy"
        self["energy_energy"].append(energy)
        self["energy_basepairing_x"] = [base_pairing_x]
        self["energy_basepairing_y"] = [base_pairing_y]

        return energy


    def convert_bracket_to_numbered_pairs(self,bracket_string):
        #all_seq_len = len(bracket_string)
        bp_x = []
        bp_y = []
        strands = []

        for y in range(bracket_string.count(")")):
            bp_y.append([])

        last_nt_x_list = []
        counter=0
        num_strands=0
        for (pos,letter) in enumerate(bracket_string[:]):
            if letter is ".":
                counter += 1

            elif letter is "(":
                bp_x.append(pos-num_strands)
                last_nt_x_list.append(pos-num_strands)
                counter += 1

            elif letter is ")":
                nt_x = last_nt_x_list.pop()
                nt_x_pos = bp_x.index(nt_x)
                bp_y[nt_x_pos] = pos-num_strands
                counter += 1

            elif letter is "&":
                strands.append(counter)
                counter=0
                num_strands+=1

            else:
                raise Exception("Error! Invalid character in bracket notation.")

        if len(last_nt_x_list) > 0:
            raise Exception("Error! Leftover unpaired nucleotides when converting from bracket notation to numbered base pairs.")

        strands.append(counter)
        bp_x = [pos+1 for pos in bp_x[:]] #Shift so that 1st position is 1
        bp_y = [pos+1 for pos in bp_y[:]] #Shift so that 1st position is 1

        return (strands,bp_x, bp_y)


    def convert_numbered_pairs_to_bracket(self,strands,bp_x,bp_y):
        bp_x = [pos-1 for pos in bp_x[:]] #Shift so that 1st position is 0
        bp_y = [pos-1 for pos in bp_y[:]] #Shift so that 1st position is 0

        bracket_notation = []
        counter=0
        for (strand_number,seq_len) in enumerate(strands):
            if strand_number > 0: bracket_notation.append("&")
            for pos in range(counter,seq_len+counter):
                if pos in bp_x:
                    bracket_notation.append("(")
                elif pos in bp_y:
                    bracket_notation.append(")")
                else:
                    bracket_notation.append(".")
            counter+=seq_len

        return "".join(bracket_notation)


    def _run_command(self, cmd, input_string):
        """Run the specified command and return output"""
        p = Popen(cmd, stdout=PIPE, stdin=PIPE, stderr=PIPE, universal_newlines=True)
        out, stderr = p.communicate(input=input_string)
        if p.returncode:
            raise Exception("Cmd {} failed: {}".format(cmd[0], stderr))
        return out
