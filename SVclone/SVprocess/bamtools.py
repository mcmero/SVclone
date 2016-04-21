'''
Copyright (c) 2013, Andreas Heger, MRC CGAT

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
    Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in
    the documentation and/or other materials provided with the
    distribution.  Neither the name of the Medical Research Council nor the
    names of its contributors may be used to endorse or promote
    products derived from this software without specific prior written
    permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

BamTools - utilities for working with BAM files
===============================================
'''

import numpy
import pysam


def isPaired(bamfile, alignments=1000):
    '''check if a *bamfile* contains paired end data

    The method reads at most the first *alignments* and returns
    True if any of the alignments are paired.
    '''

    samfile = pysam.AlignmentFile(bamfile)
    n = 0
    for read in samfile:
        if read.is_paired:
            break
        n += 1
        if n == alignments:
            break

    samfile.close()

    return n != alignments


def estimateInsertSizeDistribution(bamfile,
                                   alignments=50000):
    '''estimate insert size from first alignments in bam file.

    returns mean and stddev of insert sizes.
    '''

    assert isPaired(bamfile), \
        'can only estimate insert size from' \
        'paired bam files'

    samfile = pysam.AlignmentFile(bamfile)
    # only get positive to avoid double counting
    inserts = numpy.array(
        [read.tlen for read in samfile.head(alignments)
         if read.is_proper_pair and read.tlen > 0])
    return numpy.mean(inserts), numpy.std(inserts)


def estimateTagSize(bamfile,
                    alignments=10,
                    multiple="error"):
    '''estimate tag size from first alignments in file.'''
    samfile = pysam.AlignmentFile(bamfile)
    sizes = [read.rlen for read in samfile.head(alignments)]
    mi, ma = min(sizes), max(sizes)

    if mi == 0 and ma == 0:
        sizes = [read.inferred_length for read in samfile.head(alignments)]
        # remove 0 sizes (unaligned reads?)
        sizes = [x for x in sizes if x > 0]
        mi, ma = min(sizes), max(sizes)

    if mi != ma:
        if multiple == "error":
            raise ValueError('multiple tag sizes in %s: %s' % (bamfile, sizes))
        elif multiple == "mean":
            mi = int(sum(sizes) / len(sizes))

    return mi


def getNumberOfAlignments(bamfile):
    '''
    return number of alignments in bamfile.
    '''
    samfile = pysam.AlignmentFile(bamfile)
    return samfile.mapped


