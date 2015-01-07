import matplotlib
matplotlib.use('Agg')

import collections
import os
import sys
import re
import vcf
import numpy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from fractions import Fraction
import scipy.stats
import sklearn.linear_model

import argparse

def calc_allelic_imbalance(read_count_1, read_count_2):
    total_count = read_count_1 + read_count_2
    if read_count_1 > total_count/2.0:
        return scipy.stats.binom.sf(read_count_1, total_count, 0.5)
    else:
        return scipy.stats.binom.cdf(read_count_1, total_count, 0.5)

# Returns the index type of the locus genotype
# 0: Homozygous   reference
# 1: Homozygous   non-reference
# 2: Heterozygous reference/non-reference
# 3: Heterozygous non-reference/non-reference
def categorize_locus(a1, a2):
    if (a1 == a2):
        # Homozygous            
        if a1 == 0:
            # Reference
            return 0
        else:
            # Non-reference
            return 1
    else:
        # Heterozygous            
        if (a1 == 0) or (a2 == 0):
            # Reference/Non-reference
            return 2
        else:
            # Non-reference/Non-reference
            return 3

def get_locus_type_names():
    return ["Hom. ref", "Hom. non-ref", "Het. ref/non-ref", "Het. non-ref/non-ref"]

class CALC_OPTIMAL_OFFSETS:
    def __init__(self, marshfield_reader, output_file, min_samp_count=20):
        self.marshfield_reader = marshfield_reader
        self.output = open(output_file, "w")
        self.min_samp_count = min_samp_count
      
    def initialize(self, vcf_reader):
        sample_indices = dict(map(lambda x: reversed(x), enumerate(vcf_reader.samples)))
        self.samples, self.indexes = [], []
        for s in self.marshfield_reader.samples:
            if s in sample_indices:
                self.samples.append(s)
                self.indexes.append(sample_indices[s])
        self.comp_count = 0
        self.skip_count = 0
          
    def process_record(self, record):
        locus_id = record.CHROM.replace("chr", "") + ":" + str(record.POS)
        if not self.marshfield_reader.has_locus(locus_id):
            self.skip_count += 1
            return
        
        self.comp_count += 1
        
        # Compute the set of possible corrections
        obs_diffs  = set([])
        samp_count = 0
        for i in xrange(len(self.indexes)):
            info = record.samples[self.indexes[i]]
            if info['GT'] is None:
                continue

            marsh_gts = self.marshfield_reader.get_gts(locus_id, self.samples[i])
            other_gbs = map(int, info['GB'].split("/"))

            # Missing Marshfield data
            if marsh_gts == [0, 0]:
                continue
            obs_diffs.add(marsh_gts[0] - other_gbs[0])
            obs_diffs.add(marsh_gts[0] - other_gbs[1])
            obs_diffs.add(marsh_gts[1] - other_gbs[0])
            obs_diffs.add(marsh_gts[1] - other_gbs[1])
            samp_count += 1
        
        if samp_count < self.min_samp_count:
            return
            
        # Determine the optimal correction
        diff_scores = collections.defaultdict(int)
        for i in xrange(len(self.indexes)):
            info = record.samples[self.indexes[i]]
            if info['GT'] is None:
                continue
        
            bp_1a, bp_1b = map(int, info['GB'].split("/"))
            bp_2a, bp_2b = self.marshfield_reader.get_gts(locus_id, self.samples[i])
            
            # Missing Marshfield data
            if bp_2a == 0 and bp_2b == 0:
                continue
            
            for diff in obs_diffs:
                bp_a = bp_2a - diff
                bp_b = bp_2b - diff

                if bp_1a == bp_1b:
                    if bp_a == bp_b:
                        score = 1.0 if bp_a == bp_1a else 0.0
                    else:
                        score = 0.5 if (bp_a == bp_1a or bp_b == bp_1a) else 0.0
                else:
                    if bp_a == bp_b:
                        score = 0.25 if (bp_a == bp_1a or bp_a == bp_1b) else 0.0
                    else:
                        score = 0.0
                        if bp_a == bp_1a or bp_b == bp_1a:
                            score += 0.25
                        if bp_a == bp_1b or bp_b == bp_1b:
                            score += 0.25
                        if score == 0.5:
                            score = 1.0
                diff_scores[diff] += score
                
        # Write the optimal correction to an output file
        if len(diff_scores) > 0:
            opt_diff = sorted(diff_scores.items(), key=lambda x:x[1])[-1][0]
            self.output.write(locus_id + "\t" + str(opt_diff) + "\n")
            
    def finish(self):
        print(self.comp_count, self.skip_count)
        self.output.close()
    

class CALC_OPTIMAL_OFFSETS_OTHER:
    def __init__(self, marshfield_reader, output_file):
        self.marshfield_reader = marshfield_reader
        self.output = open(output_file, "w")
      
    def initialize(self, vcf_reader):
        sample_indices = dict(map(lambda x: reversed(x), enumerate(vcf_reader.samples)))
        self.samples, self.indexes = [], []
        for s in self.marshfield_reader.samples:
            if s in sample_indices:
                self.samples.append(s)
                self.indexes.append(sample_indices[s])
        self.comp_count, self.skip_count = 0, 0
        
    def process_record(self, record):
        locus_id = record.CHROM.replace("chr", "") + ":" + str(record.POS)
        if not self.marshfield_reader.has_locus(locus_id):
            self.skip_count += 1
            return
        
        self.comp_count += 1
            
        # Determine the optimal correction
        diff_scores = collections.defaultdict(int)
        for i in xrange(len(self.indexes)):
            info = record.samples[self.indexes[i]]
            if info['GT'] is None:
                continue
            
            bp_1a, bp_1b = map(int, info['GB'].split("/"))
            bp_2a, bp_2b = self.marshfield_reader.get_gts(locus_id, self.samples[i])
            
            # Missing Marshfield data
            if bp_2a == 0 and bp_2b == 0:
                continue
            
            if bp_2a == bp_2b:
                if bp_1a == bp_1b:
                    diff_scores[bp_2a-bp_1a] += 1.0
                else:
                    diff_scores[bp_2a-bp_1a] += 0.5
                    diff_scores[bp_2a-bp_1b] += 0.5
            else:
                if bp_1a == bp_1b:
                    diff_scores[bp_2a-bp_1a] += 0.5
                    diff_scores[bp_2b-bp_1a] += 0.5
                else:
                    d11 = bp_2a - bp_1a
                    d12 = bp_2b - bp_1a
                    d21 = bp_2a - bp_1b
                    d22 = bp_2b - bp_1b
                    
                    if d11 == d22:
                        diff_scores[d11] += 1.0
                    elif d12 == d21:
                        diff_scores[d12] += 1.0
                    else:
                        diff_scores[d11] += 0.5
                        diff_scores[d12] += 0.5
            
        if len(diff_scores) > 0:
            opt_diff = sorted(diff_scores.items(), key=lambda x:x[1])[-1][0]
            self.output.write(locus_id + "\t" + str(opt_diff) + "\n")
            
    def finish(self):
        print(self.comp_count, self.skip_count)
        self.output.close()


class MARSHFIELD_COMPARATOR:
    def __init__(self, marshfield_reader, corrections={}, norm_by_period=False, scale_by_width=False, plot_line=True, inc_ref_len=False, debug_file=None):
        self.marshfield_reader = marshfield_reader
        self.corrections       = corrections
        self.norm_by_period    = norm_by_period
        self.scale_by_width    = scale_by_width
        self.plot_line         = plot_line
        self.debug_file        = debug_file
        self.inc_ref_len       = inc_ref_len
        
    def initialize(self, vcf_reader):
        sample_indices = dict(map(lambda x: reversed(x), enumerate(vcf_reader.samples)))
        self.samples, self.indexes = [], []
        for s in self.marshfield_reader.samples:
            if s in sample_indices:
                self.samples.append(s)
                self.indexes.append(sample_indices[s])
        
        self.total_type_counts    = numpy.zeros(4)
        self.total_correct_counts = numpy.zeros(4)
        self.pair_counts          = collections.defaultdict(int)
        self.x_vals, self.y_vals  = [], []
        self.proc_samples         = set([])
        self.num_proc_calls       = 0        
        self.call_types           = [[0, 0, 0, 0],       [0, 0, 0, 0],      [0, 0, 0, 0, 0],         [0, 0, 0, 0, 0]]
        self.lobstr_vals          = [[[], [], [], []],   [[], [], [], []],   [[], [], [], [], []],   [[], [], [], [], []]]
        self.marsh_vals           = [[[], [], [], []],   [[], [], [], []],   [[], [], [], [], []],   [[], [], [], [], []]]
        self.lobstr_vals_diff     = [[[], [], [], []],   [[], [], [], []],   [[], [], [], [], []],   [[], [], [], [], []]]
        self.marsh_vals_diff      = [[[], [], [], []],   [[], [], [], []],   [[], [], [], [], []],   [[], [], [], [], []]]
        self.debug_locus_info     = {}
        self.accuracy_per_locus   = {}
        self.correct_indicators   = []
        self.read_depths          = []
        self.quality_scores       = []
        self.heterozygotes        = []

    def process_record(self, record):
        locus_id = record.CHROM.replace("chr", "") + ":" + str(record.POS)
        if not self.marshfield_reader.has_locus(locus_id):
            return
        if locus_id in self.corrections:
            correction = self.corrections[locus_id]
        else:
            return
        
        locus_process_count, locus_correct_count = 0, 0
        for i in xrange(len(self.indexes)):
            info = record.samples[self.indexes[i]]
            if info['GT'] is None:
                continue            
            if self.marshfield_reader.get_gts(locus_id, self.samples[i]) == [0, 0]:
                continue
                
            locus_process_count += 1
            self.proc_samples.add(self.samples[i])
            self.num_proc_calls += 1
            
            gb_1a, gb_1b = map(lambda x: x-correction, self.marshfield_reader.get_gts(locus_id, self.samples[i]))
            gb_2a, gb_2b = map(int, info['GB'].split("/"))
            total_one    = gb_1a + gb_1b
            total_two    = gb_2a + gb_2b
            if self.inc_ref_len:
                total_one += 2*len(record.alleles[0])
                total_two += 2*len(record.alleles[0])

            locus_type = categorize_locus(gb_1a, gb_1b)
            self.total_type_counts[locus_type] += 1
            
            if locus_type == 0 or locus_type == 1:
                # Homzygous ref or homzygous non-ref allele
                type_index = 0 if gb_1a == 0 else 1
                    
                if gb_2a == gb_2b:
                    sub_index = 0 if gb_2a == gb_1a else 1
                else:
                    sub_index = 2 if (gb_2a == gb_1a or gb_2b == gb_1a) else 3
            else:
                # Heterozygous 
                type_index = 2 if (gb_1a == 0 or gb_1b == 0) else 3
                    
                if gb_2a == gb_2b:
                    sub_index = 0 if (gb_2a == gb_1a or gb_2a == gb_1b) else 1
                else:
                    corr_count = 0 + (gb_2a == gb_1a) + (gb_2a == gb_1b) + (gb_2b == gb_1a) + (gb_2b == gb_1b)
                    sub_index  = corr_count + 2

            if (gb_1a == gb_2a and gb_1b == gb_2b) or (gb_1a == gb_2b and gb_1b == gb_2a):
                locus_correct_count += 1
                self.correct_indicators.append(1)
            else:
                self.correct_indicators.append(0)
            self.read_depths.append(info['DP'])
            self.quality_scores.append(info['Q'])
            self.heterozygotes.append(gb_1a != gb_1b)

            should_record = False
            if locus_type == 0 or locus_type == 1:
                if sub_index == 2 or sub_index == 3:        
                    should_record = True
            else:
                if sub_index == 1 or sub_index == 2 or sub_index == 3:
                    should_record = True
            if should_record:
                chrom, start = locus_id.split(":")
                key = ("chr"+chrom, start, str(int(start)+len(record.alleles[0])))
                if key not in self.debug_locus_info:
                    self.debug_locus_info[key] = []
                self.debug_locus_info[key].append([self.samples[i], (gb_1a, gb_1b), (gb_2a, gb_2b)])

            self.call_types[type_index][sub_index] += 1
            self.marsh_vals[type_index][sub_index].append( gb_1a+gb_1b+2*len(record.alleles[0]))
            self.lobstr_vals[type_index][sub_index].append(gb_2a+gb_2b+2*len(record.alleles[0]))
            self.marsh_vals_diff[type_index][sub_index].append( gb_1a+gb_1b)
            self.lobstr_vals_diff[type_index][sub_index].append(gb_2a+gb_2b)
            
            if total_one == total_two:
               self.total_correct_counts[locus_type] += 1
               
            if self.norm_by_period:
                motif = record.INFO['MOTIF']
                key = (Fraction(total_one, len(motif)), Fraction(total_two, len(motif)))
                self.x_vals.append(1.0*total_one/len(motif))
                self.y_vals.append(1.0*total_two/len(motif))
            else:
                key = (total_one, total_two)
                self.x_vals.append(total_one)
                self.y_vals.append(total_two)
            self.pair_counts[key] += 1

            # Print line for database
            # print("%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s"%(gb_2a, gb_2b, gb_1a, gb_1b, total_one, total_two, locus_id.split(":")[0], locus_id.split(":")[1], self.samples[i]))

        if locus_process_count != 0:
            self.accuracy_per_locus[locus_id] = [locus_correct_count, locus_process_count, 1.0*locus_correct_count/locus_process_count]

            
    def finish(self, pdfpage):
        for pair in sorted(self.accuracy_per_locus.items(), key = lambda x: x[1][2]):
            print(pair[0], pair[1])

        for valset in [self.read_depths, self.quality_scores]:
            pairs = sorted(zip(valset, self.correct_indicators))
            x, y, z    = [], [], []
            prev_val   = pairs[-1][0]+1
            correct    = 0
            total      = 0
            for i in xrange(len(pairs)-1, -1, -1):
                correct += pairs[i][1]
                total   += 1
                if pairs[i][0] != prev_val:
                    x.append(pairs[i][0])
                    y.append(100.0*correct/total)
                    z.append(total)
            fig = plt.figure()
            ax  = fig.add_subplot(111)
            ax.plot(x, y)
            ax.set_xlabel("Threshold")
            ax.set_ylabel("% Correct")
            ax2 = ax.twinx()
            ax2.plot(x, z, 'r')
            pdfpage.savefig(fig)
        
        items    = filter(lambda x: x[2], sorted(zip(self.read_depths, self.correct_indicators, self.heterozygotes)))
        index    = 0
        depths   = []
        num_corr = []
        num_tot  = []
        while index < len(items):
            prev_depth  = items[index][0]
            num_correct = 0
            num_total   = 0
            while index < len(items) and items[index][0] == prev_depth:
                num_correct += items[index][1]
                num_total   += 1
                index       += 1
            depths.append(prev_depth)
            num_corr.append(num_correct)
            num_tot.append(num_total)

        print(depths, num_corr, num_tot)

        fig = plt.figure()
        ax  = fig.add_subplot(111)
        ax.plot(depths, 1.0*numpy.array(num_corr)/numpy.array(num_tot), 'g-o')
        ax.set_xlabel("Minimum Read Depth")
        ax.set_ylabel("% Heterozygotes Correct")
        ax2 = ax.twinx()
        ax2.plot(depths, num_tot, 'r-o')
        ax2.set_yscale('log')
        pdfpage.savefig(fig)

        fig = plt.figure()
        ax  = fig.add_subplot(111)
        num_corr = numpy.cumsum(num_corr[::-1])[::-1]
        num_tot  = numpy.cumsum(num_tot[::-1])[::-1]
        print(depths, num_corr, num_tot)

        ax.plot(depths, 1.0*numpy.array(num_corr)/numpy.array(num_tot), 'g-o')
        ax.set_xlabel("Minimum Read Depth")
        ax.set_ylabel("% Heterozygotes Correct")
        ax2 = ax.twinx()
        ax2.plot(depths, num_tot, 'r-o')
        ax2.set_yscale('log')
        pdfpage.savefig(fig)

        fig, axes = plt.subplots(2, 1, sharex=True, sharey=False)
        axes[0].plot(depths, 100.0*numpy.array(num_corr)/numpy.array(num_tot), '-o')
        axes[0].set_ylabel("% Heterozygotes correct")
        axes[1].plot(depths, num_tot, '-o')
        axes[1].set_yscale('log')
        axes[1].set_xlabel("Minimum read depth")
        axes[1].set_ylabel("Number of relevant calls")
        axes[0].xaxis.set_ticks_position('bottom')
        axes[0].yaxis.set_ticks_position('left')
        axes[0].set_ylim((0, 105))
        axes[1].xaxis.set_ticks_position('bottom')
        axes[1].yaxis.set_ticks_position('left')
        pdfpage.savefig(fig)
      

        if self.debug_file is not None:
            debug_stream = open(self.debug_file, "w")
            for locus,vals in self.debug_locus_info.items():
                samples = map(lambda x: x[0], vals)
                annots  = map(lambda x: str(x[1])+"|"+str(x[2]), vals)
                annots  = map(lambda x: x.replace(", ","/"), annots)
                sample_string = ",".join(samples)
                annot_string  = ",".join(annots)
                debug_stream.write("%s\t%s\t%s\t%s\t%s\n"%(locus[0], locus[1], locus[2], sample_string, annot_string))
            debug_stream.close()

        print("Processed a total of %d loci"%len(self.accuracy_per_locus))
        print("Processed a total of %d samples"%len(self.proc_samples))
        print("Processed a total of %d calls"%self.num_proc_calls)
        
        # Compute concordance r^2
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(self.x_vals, self.y_vals)

        new_x = numpy.vstack(numpy.array(self.x_vals))
        new_y = numpy.array(self.y_vals).T
        print(new_x.shape, new_y.shape)

        L1_model = sklearn.linear_model.ElasticNet(alpha=1.0, l1_ratio=1.0)
        L1_model.fit(new_x, new_y)
        print(L1_model.coef_, L1_model.intercept_)

        print("Best-fit line stats: y = %f*x + %f"%(slope, intercept))
        print("Concordance R^2= "+str(r_value*r_value))
        print("Correctness Stats:")
        print("\t                (Hom.Ref, Hom. Non-ref, Het. w/Ref, Het w/o Ref)")
        print("\tTotal  Correct: "+str(self.total_correct_counts))
        print("\tTotal  Analyzed:"+str(self.total_type_counts))

        A     = numpy.vstack([self.x_vals]).T
        slope = numpy.linalg.lstsq(A, self.y_vals)[0][0]
        print("\tSlope = %f when forcing a zero intercept"%(slope))

        # Print concordance statistics by type
        frac_correct = 1.0*self.total_correct_counts/self.total_type_counts
        type_names   = get_locus_type_names()
        print("Correctness statistics by type:")
        for i in xrange(4):
            print("\t"+type_names[i]+"\t"+str(frac_correct[i]))
        print("\tOverall \t" + str(100.0*numpy.sum(self.total_correct_counts)/numpy.sum(self.total_type_counts)))

        print("Call type stats:")
        type_names = ["Hom. Ref", "Hom. Non-Ref", "Het. w/Ref", "Het. w/o Ref"]
        for i in xrange(4):
            print(type_names[i] + "\t" + str(self.call_types[i]) + "\t" + str(1.0*numpy.array(self.call_types[i])/sum(self.call_types[i])))
        print("\n")

        print("R^2 stats:")
        for i in xrange(4):
            rsq_vals = []
            for j in xrange(len(self.lobstr_vals[i])):
                if len(self.lobstr_vals[i][j]) > 2:
                    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(self.lobstr_vals[i][j], self.marsh_vals[i][j])
                    rsq_vals.append(r_value*r_value)
                else:
                    rsq_vals.append("NA")
            print(type_names[i] + "\t" + str(rsq_vals))
            
        print("R^2 for total allele lengths:")
        for i in xrange(4):
            x = reduce(lambda x,y: x+y, self.lobstr_vals[i])
            y = reduce(lambda x,y: x+y, self.marsh_vals[i])
            if len(x) > 2:
                slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
                print(type_names[i] + "\t" + str(r_value*r_value))
            else:
                print(type_names[i] + "\t" + "NA")

        print("R^2 for difference from reference:")
        for i in xrange(4):
            x = reduce(lambda x,y: x+y, self.lobstr_vals_diff[i])
            y = reduce(lambda x,y: x+y, self.marsh_vals_diff[i])
            if len(x) > 2:
                slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
                print(type_names[i] + "\t" + str(r_value*r_value))
            else:
                print(type_names[i] + "\t" + "NA")

        homozygous_correct   = numpy.array(self.call_types[0]) + numpy.array(self.call_types[1])
        homozygous_total     = sum(self.call_types[0]) + sum(self.call_types[1])
        heterozygous_correct = numpy.array(self.call_types[2]) + numpy.array(self.call_types[3])
        heterozygous_total   = sum(self.call_types[2]) + sum(self.call_types[3])

        # Arrange data for bubble plot
        x    = []
        y    = []
        vals = []
        for k, v in self.pair_counts.iteritems():
            x.append(float(k[0]))
            y.append(float(k[1]))
            vals.append(1.0*v)

        # Create bubble plot
        fig = plt.figure()
        ax  = fig.add_subplot(111)
        
        if self.inc_ref_len:
            ax.set_xlabel("Sum of Marshfield allele lengths")
            ax.set_ylabel("Sum of lobSTR allele lengths") 
        else:
            ax.set_xlabel("Marshfield Genotype Sum")
            ax.set_ylabel("lobSTR Genotype Sum")
            
        # color=(1.0*41/255, 1.0*171/255, 1.0*226/255)
        if self.scale_by_width:
            ax.scatter(x, y, s=numpy.pi*numpy.array(vals)*numpy.array(vals)/4.0, alpha=0.59)
        else:
            ax.scatter(x, y, s=numpy.array(vals), alpha=0.59)
        
        ax.set_alpha(0.8)
        min_lim = min(min(x), min(y))
        max_lim = max(max(x), max(y))
        lim     = max(-min_lim, max_lim)

        lim = 80
        #lim = 40

        if self.plot_line:
            if self.inc_ref_len:
                line_x = numpy.arange(0, 150, 1.0)
            else:
                line_x = numpy.arange(-lim+10, lim-10, 1.0)
            line_y = line_x
            ax.plot(line_x, line_y, linestyle='--', color='k')

        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        if self.inc_ref_len:
            ax.set_xlim(0, 150)
            ax.set_ylim(0, 150)
            
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.set_aspect('equal')
        pdfpage.savefig(fig)

        fig = plt.figure()
        ax  = fig.add_subplot(111)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.set_aspect('equal')
        leg_x    = [0,  0,  0,   0]
        leg_y    = [0, lim/4.0, lim/2.0, 3*lim/4.0]
        #text_x   = [0,  0,  0,   0]
        text_x   = [lim/4.0, lim/4.0, lim/4.0, lim/4.0]

        text_y   = [0, lim/4.0, lim/2.0, 3*lim/4.0]
        #text_y   = [lim/8.0, 3*lim/8.0, 5*lim/8.0, 7*lim/8.0]
        leg_size = [25, 100, 500]
        for x,y,size,tx,ty in zip(leg_x, leg_y, leg_size, text_x, text_y):
            if self.scale_by_width:
                ax.scatter([x], [y], s=numpy.pi*numpy.array(size)*numpy.array(size)/4.0) 
            else:
                ax.scatter([x], [y], s=numpy.array([size]))
            ax.text(tx, ty, str(size), {'ha':'center', 'va':'center'})
        pdfpage.savefig(fig)
        
        '''
        print("Individual datapoints:")
        for i in xrange(len(self.x_vals)):
            print("%d\t%d"%(self.x_vals[i], self.y_vals[i]))
        '''
            
class MARSHFIELD_READER:
    def __init__(self, input_file):
        data         = open(input_file, "r")
        self.samples = data.readline().strip().split("\t")[10:]
        self.sample_indices = {}
        for i,s in enumerate(self.samples):
            self.sample_indices[s] = i
        
        self.locus_indices = {}
        self.gt_data       = []
        locus_count        = 0
        for line in data:
            tokens = line.strip().split("\t")
            key    = "%s:%s"%(tokens[0].replace("chr",""), tokens[1])
            gts    = map(lambda x: map(int, x.split("/")), tokens[10:])
            self.gt_data.append(gts)
            self.locus_indices[key] = locus_count
            locus_count += 1
        data.close()
        
    def has_locus(self, key):
        return key in self.locus_indices

    def has_sample(self, sample):
        return sample in self.sample_indices
        
    def get_gts(self, locus, sample):
        return self.gt_data[self.locus_indices[locus]][self.sample_indices[sample]]
    
class AVG_COVERAGE:
    def initialize(self, dbase):
        self.coverage_index = dbase.format_fields.get_field_index("DP")
    def process_line(self, line, dbase):
        sum_coverage = 0
        num_samples  = 0
        for samp_index in dbase.sample_indexes.values():
            if dbase.missing_data[samp_index]:
                continue
            num_samples  += 1
            sum_coverage += dbase.get_sample_data(samp_index, self.coverage_index)
        if num_samples == 0:
            return 0
        else:
            return 1.0*sum_coverage/num_samples
    def finish(self, pdfpage):
        pass

class FRAC_MISSING:
    def initialize(self, dbase):
        pass
    def process_line(self, line, dbase):
        missing_count = 0
        for samp_index in dbase.sample_indexes.values():
            if dbase.missing_data[samp_index]:
                missing_count += 1
        return 1.0*missing_count/len(dbase.sample_indexes)
    def finish(self, pdfpage):
        pass

class REF_LENGTH:
    def initialize(self, dbase):
        pass
    def process_line(self, line, dbase):
        return len(dbase.reference)
    def finish(self, pdfpage):
        pass

class LOCUS_FILTER:
    def __init__(self, processor, filter_fxn):
        self.processor  = processor
        self.filter_fxn = filter_fxn
    def initialize(self, dbase):
        self.processor.initialize(dbase)
    def process_line(self, line, dbase):
        return self.filter_fxn(self.processor.process_line(line, dbase))
    def finish(self, pdfpage):
        pass
    
def read_corrections(input_file):
    corrections = {}
    data = open(input_file, "r")
    for line in data:
        tokens = line.strip().split()
        corrections[tokens[0]] = int(tokens[1])
    data.close()
    return corrections

def main():
    filter_low_counts = False
    filter_dropouts   = False
    filter_imbalanced = False

    parser = argparse.ArgumentParser()
    parser.add_argument("--corr",     required=True,  help="Path for BED file of corrections")
    parser.add_argument("--calls",    required=True,  help="File containing calls to compare")
    parser.add_argument("--marsh",    required=True,  help="Marshfield call file")
    parser.add_argument("--out",      required=True,  help="Output file")
    parser.add_argument("--debug",    required=False, default=None,                help="")
    
    '''
    parser.add_argument("--min_miss", required=False, type=float, default=0.0,     help="Minimum fraction of missing calls")
    parser.add_argument("--max_miss", required=False, type=float, default=1.0,     help="Maximum fraction of missing calls")
    parser.add_argument("--min_cov",  required=False, type=float, default=0.0,     help="Minimum average coverage")
    parser.add_argument("--max_cov",  required=False, type=float, default=10000.0, help="Maximum average coverage")
    parser.add_argument("--min_len",  required=False, type=int,   default=0,       help="Minimum reference length")
    parser.add_argument("--max_len",  required=False, type=int,   default=10000,   help="Maximum reference length")
    '''
    parser.add_argument("--min_samp_for_corr", required=False, type=int, default=20, help="Minimum number of samples required to estimate a length correction")
    args = parser.parse_args()
    
    marsh_file   = args.marsh
    lobSTR_file  = args.calls
    output_file  = args.out
    corr_file    = args.corr
    
    # Calculate the optimal corrections
    marsh_reader = MARSHFIELD_READER(marsh_file)
    
    if not os.path.exists(corr_file):
        print("Recalculating optimal offsets")
        offset_calc = CALC_OPTIMAL_OFFSETS(marsh_reader, corr_file, args.min_samp_for_corr)
        vcf_reader  = vcf.Reader(filename=lobSTR_file)
        offset_calc.initialize(vcf_reader)
        for record in vcf_reader:
            offset_calc.process_record(record)
        offset_calc.finish()
    else:
        print("Using preexisting offsets")        

    # Assess the concordance and create a bubble plot
    corrections = read_corrections(corr_file)
    marsh_comp  = MARSHFIELD_COMPARATOR(marsh_reader, corrections=corrections, norm_by_period=False, scale_by_width=False, plot_line=True, inc_ref_len=True, debug_file=args.debug)

    #dbase.add_filter(LOCUS_FILTER(REF_LENGTH(),   lambda x: x >= args.min_len))
    #dbase.add_filter(LOCUS_FILTER(REF_LENGTH(),   lambda x: x <= args.max_len))
    #dbase.add_filter(LOCUS_FILTER(FRAC_MISSING(), lambda x: x >= args.min_miss))
    #dbase.add_filter(LOCUS_FILTER(FRAC_MISSING(), lambda x: x <= args.max_miss))
    #dbase.add_filter(LOCUS_FILTER(AVG_COVERAGE(), lambda x: x >= args.min_cov))
    #dbase.add_filter(LOCUS_FILTER(AVG_COVERAGE(), lambda x: x <= args.max_cov)) 

    pp = PdfPages(output_file)
    vcf_reader = vcf.Reader(filename=lobSTR_file)
    marsh_comp.initialize(vcf_reader)
    for record in vcf_reader:
        marsh_comp.process_record(record)
    marsh_comp.finish(pp)
    pp.close()
    return



if __name__ == "__main__":
    main()
