import os
import csv
import unittest
import re
import numpy as np
import sv_proc
from unittest import TestCase

class test_proc_svs(unittest.TestCase):
    # do stuff...
    def test_is_same_sv(self):
        self.assertTrue(phyl_sv.is_same_sv(x[0],y[0],ff))
        self.assertTrue(phyl_sv.is_same_sv(x[0],y[1],ff))
        self.assertFalse(phyl_sv.is_same_sv(x[2],y[0],ff))
        self.assertTrue(phyl_sv.is_same_sv(x[3],y[3],ff))
