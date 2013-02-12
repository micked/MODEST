#!/usr/bin/env python

from mage_tool.helpers import o2n
from mage_tool.IO import Mutation

if __name__ == "__main__":
    print(o2n(5, "AA", "Insertion", "AAAACCGG"))
    print(Mutation("arrow", "AA", 5, "Insertion"))
