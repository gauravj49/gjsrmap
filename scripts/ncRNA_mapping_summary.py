#!/usr/local/bin/python
"""
***********************************************
- PROGRAM: ncRNA_mapping_summary.py
- CONTACT: Gaurav Jain(gaurav.jain@dzne.edu)
***********************************************
"""
print (__doc__)

# Built in modules
import argparse
import os.path
import sys

# 3rd party modules
import textwrap
import re
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib as mp
#mp.use('Agg') # to use matplotlib without X11
mp.use('pdf') # to use matplotlib without X11
import matplotlib.pyplot as plt
from matplotlib import rcParams
import subprocess
import binascii as bi
import scipy.stats as stats
from collections import *
from numpy import nanmean

# for looping files in a dir
import glob

# user defined modules
from gjainPyLib import *      # import all the functions from the Gaurav`s python scripts/library

### for color scale
from  matplotlib import colors
from itertools import cycle, islice, chain # barplot colors
from termcolor import colored

### Modeling related ###
from time import time
from operator import itemgetter


################ USER CONFIGURATION ###################
# 1)  Save matplotlib figure with text as vectorized text and not as a path
#     The matplotlib documentation in http://matplotlib.org/users/customizing.html 
#     states that the default value for the parameter svg.fonttype is 'path', 
#     which means that characters will be converted to paths when exporting to svg format.
matplotlib.rcParams['svg.fonttype']    = 'none'
matplotlib.rcParams['pdf.fonttype']    =  42
matplotlib.rcParams['patch.edgecolor'] = 'white'
#######################################################

def main():
    # Get input options
    args               = check_options()
    input_counts_dir   = args.input_counts_dir.rstrip('\/')
    output_file        = args.output_file
    subset_files       = args.subset_files
    merged_excel       = args.merged_excel
    skip_subset_conf   = args.skip_subset_conf
    subset_sample_file = args.subset_sample_file
    min_library_size   = args.min_library_size

    if not output_file:
        output_file = "{0}/summary_plots/all/all.txt".format(input_counts_dir)
    create_dir("{0}".format(get_file_info(output_file)[0]))

    # Get the list of allncrnaCount files
    allncnorm_counts_dir = "{0}/normalized_counts/allncrna/".format(input_counts_dir) 
    allncrnaCount_files = [ f for f in os.listdir(allncnorm_counts_dir) if os.path.isfile(os.path.join(allncnorm_counts_dir, f))]

    # Create separate directories for pdfs
    output_pdf_dir = "{0}/pdf".format(get_file_info(output_file)[0])
    create_dir(output_pdf_dir)
    
    # Subset the list of files
    if subset_files or subset_sample_file:
        if subset_sample_file:
            subset_files = '|'.join([line.strip() for line in open(subset_sample_file, 'r')])
            # if the last character is '|' ex: 'A|B|' ... we want 'A|B'
            subset_files = subset_files.rstrip('|')
            subset_files = subset_files.replace("_allncrnaCounts.txt", "")

        # Create the regular expression
        r = re.compile(subset_files)
        allncrnaCount_files = list(filter(r.search, allncrnaCount_files))
        print("- Subset filenames:")
        for f in allncrnaCount_files:
            print("\t- {0}".format(f))
        
        print("\n- Total subset samples: {0}".format(len(allncrnaCount_files)))
        if not skip_subset_conf:
            # Confirm the file names from user
            print("\n- Do you want to proceed? Please press 'y' for yes or 'n' for no.")
            print("\t- Your option: ", end=' ')
            option = input()
            if option.lower() == 'n':
                sys.exit('You pressed "NO". Exiting now!')

    # Define smallrna class and its colors
    # Note: allncrnaCounts is for library size (please see get_counts() function for more details)
    smallrna_type = ['mirna', 'pirna', 'rRNA'  , 'snRNA'  , 'snoRNA' , 'premirna', 'osncRNA', 'allncrna']
    smallrna_cols = ['Blues', 'Reds' , 'Greens', 'Purples', 'Oranges', 'copper'  , 'summer_r', 'PuRd']
    smallrna_bcol = ['dodgerblue','firebrick' , 'forestgreen','darkorange',  'darkorchid', 'fuchsia', 'goldenrod']

    # m is a dictionary with:
    # 	primary_key    = ['count_file_id']
    # 	secondary_keys = ['mirna', 'piwirna', 'rrna', 'snorna', 'snrna', 'uniquely_mapped_reads']
    m = defaultdict(dict)

    # Uniquely mapped reads
    uniquely_mapped_reads_file = "{0}_uniquely_mapped_reads.txt".format(get_file_info(output_file)[3])
    if os.path.isfile(uniquely_mapped_reads_file):
        os.system("rm {0}".format(uniquely_mapped_reads_file))

    umrfo = open(uniquely_mapped_reads_file, 'a')

    # Get the counts from the counts file
    print("- Getting the Dictionary of counts for all smallRNAs...")
    for st in smallrna_type:
        if st == 'allncrna':
            #print colored("\n- Getting library size...", "green")
            print("\n- Getting library size...")
        else:
            #print colored("\n- Getting counts for: {0} ..".format(st), "green")
            print("\n- Getting counts for: {0} ..".format(st))

        for f in allncrnaCount_files:
            get_counts(input_counts_dir, f, st, m, umrfo)

    umrfo.close()
    os.system("sort -nk2,2 {0} -o {0}".format(uniquely_mapped_reads_file))

    # Remove list element "allncrna" from the list as we don't need to plot it
    try:
        smallrna_type.remove('allncrna')
    except Exception as e: 
        print(str(e))
        pass

    # Convert nested dictionary into pandas dataframe
    alzcltdata = pd.DataFrame.from_dict(m).transpose()

    # Dictionary keys:
    # mirna
    # mirna_total_counts
    # pirna
    # pirna_total_counts
    # osncRNA
    # osncRNA_total_counts
    # rRNA
    # rRNA_total_counts
    # snRNA
    # snRNA_total_counts
    # snoRNA
    # snoRNA_total_counts
    # uniquely_mapped_reads

    # Get the merged counts file
    print("\n- Getting merged counts file for all samples...")
    merged_file_dir = "{0}/merged_counts".format(get_file_info(output_file)[0])
    create_dir(merged_file_dir)

    if merged_excel:
        # Initialize the excel object
        merged_counts_excel = "{0}/{1}_merged_counts_normalized.xlsx".format(merged_file_dir, get_file_info(output_file)[1])
        writer = pd.ExcelWriter(merged_counts_excel)
 
    for smtype in smallrna_type:
        try:
            merged_df = pd.DataFrame()
            for f in allncrnaCount_files:
                # Get the fileid
                fileid = re.sub('_allncrnaCounts_normalized.txt', '', f)

                smf = "{0}/normalized_counts/{1}/{2}_{1}Counts_normalized.txt".format(input_counts_dir, smtype, fileid) 
                if os.path.isfile(smf):
                    # Get pandas dataframe for each fileid for each smtype
                    dfm = pd.DataFrame.from_dict(alzcltdata[smtype][fileid], orient='index')
    
                    # Add column name
                    dfm.columns = [fileid]

                    # Concatenate to the existing dataframe
                    merged_df = pd.concat([merged_df, dfm], axis=1)
                else:
                    print(colored("\t- SampleId: {0} for {1} \n\t\t- Counts file does not exists...".format(fileid, smtype), "red"))
        except Exception as e: 
            print(str(e))
            #pass

        # Add the mean as the last column
        merged_df = pd.concat([merged_df, pd.DataFrame(merged_df.mean(axis=1), columns=['mean'])], axis=1)

        # Rearrange the columns so that means is second column after featureid
        merged_cols = merged_df.columns.tolist()          # Get columns list
        merged_cols = merged_cols[-1:] + merged_cols[:-1] # move last column to first column
        merged_df   = merged_df[merged_cols]              # Reorder the dataframe with new order
        
        # Save the merged counts file
        merged_counts_file = "{0}/{1}_{2}MergedCounts_normalized.txt".format(merged_file_dir, get_file_info(output_file)[1], smtype)
        print("\t- {0}".format(merged_counts_file))
        
        # Save the text file
        fo = open(merged_counts_file, 'w')
        merged_df.to_csv(fo, sep="\t", index=True, index_label="featureid")
        fo.close()

        # Save the excel file
        if merged_excel:
            merged_df.to_excel(writer, smtype)
    

    # Save the excel writer object
    if merged_excel:
        writer.save()
        print("- {0}".format(merged_counts_excel))
    
    print("\n------------------------------------ PLOTS ------------------------------------------\n")
    print("1) Plotting library size for the samples ...")
    df = pd.Series.to_frame(alzcltdata['uniquely_mapped_reads'])
    df = df.sort_values(by = ['uniquely_mapped_reads'], ascending=[True])
    plot_library_sizes(df, output_file, min_library_size, output_pdf_dir)
    
    #-----------------------------------------------------------------------------------------------------------------------------#
    print("\n2) Plotting the aggregate of all the smallRNA classes in a pie chart...")
    # Get the aggregate of all the small rna classes
    allsmallrna = []
    for stype in smallrna_type:
        key = stype + '_total_counts'
        if key in alzcltdata:
            allsmallrna.append(np.sum(alzcltdata[key]))
        else:
            allsmallrna.append(0)
            
    # Convert the list into percentages and them make a dataframe
    smallrna_summary_df = pd.DataFrame([i*100.0/sum(allsmallrna) for i in allsmallrna], index=smallrna_type, columns=["allsmallRNA"])

    try:
        fig = plt.figure(figsize=(5,11))
        x = smallrna_summary_df.index.tolist()
        y = np.array(smallrna_summary_df["allsmallRNA"])
        pc = 100.*y/y.sum()
        patches, texts = plt.pie(y, colors=smallrna_bcol, startangle=90, radius=1.2)
        labels = ['{0:10} - {1:2.2f} %'.format(i,j) for i,j in zip(x, pc)]

        # Set the plot axis as equal to plot the pie in a circle
        plt.axis('equal')

        # Plot the legend with percentages
        plt.legend(patches, labels, loc='lower center', fontsize=12)

        plt.suptitle("Percentage of mapped SmallRNA class reads\n\n{0}".format(get_file_info(output_file)[1]))
        output_pie_plot = "{0}/{1}_piechart_smallRNA_classes.pdf".format(output_pdf_dir, get_file_info(output_file)[1])
        plt.savefig(output_pie_plot, format="pdf", bbox_inches='tight')
        output_pie_plot = get_file_info(output_file)[3] + "_piechart_smallRNA_classes.png"
        plt.savefig(output_pie_plot, format="png", bbox_inches='tight')
        
        # clear up the plot
        plt.close('all')
    except Exception as e: 
        print(str(e))
        pass

    #-----------------------------------------------------------------------------------------------------------------------------#
    # Plot most abundant smallrnas
    print("\n3) Plotting most abundant smallRNAs for every smallRNA class...")
    smallrna_abundance_barcharts(smallrna_type, smallrna_bcol, alzcltdata, 'All_nonCoding_RNA', output_file, 15)
    
    # clear up the plot
    plt.close('all')


################ USER DEFINED FUNCTIONS ###################
def plot_library_sizes(df, output_file, min_library_size, output_pdf_dir):
    ''' plot library sizes '''

    # Convert string to int
    df = df.apply(pd.to_numeric)
    df.sort_values(by=['uniquely_mapped_reads'],ascending=False, inplace=True)
    #print df.info()

    
    # Get specific colors for library size restrictions
    library_colors = []
    for i, lv in df.iterrows(): # keys are the names of the boys
        if lv['uniquely_mapped_reads'] < min_library_size:
            library_colors.append('firebrick')
        else:
            library_colors.append('black')

    nsamples  = df.shape[0]
    font_size = 20 
    if nsamples > 15:
        font_size = 10

    # Plot the bar plot
    fig = plt.figure(figsize=(nsamples, nsamples), frameon=False)
    bdf = plt.barh(list(range(len(df['uniquely_mapped_reads']))), df['uniquely_mapped_reads'], color=library_colors, align='center')

    # Get the title and ytick labels
    plt.title('Library Size\n{0}'.format(get_file_info(output_file)[1]), fontsize=font_size)
    plt.yticks(range(len(df)), df.index.tolist(), size=font_size)
    
    # Get the max and min yaxis values and let matplotlib take the log for us
    xmin, xmax = plt.gca().get_xlim()
    plt.xlim([10, xmax*1.5])
    plt.gca().set_xscale('log', basex=10)
    plt.xlabel('Log10(Uniquely mapped reads)')
    for rect in bdf:
        plt.text(1.05*rect.get_width(), rect.get_y()+0.5*rect.get_height(), '%d' % int(rect.get_width()), va='center', ha='left')

    # Create the custom legend and put below current axis
    import matplotlib.patches as mpatches
    gt = mpatches.Patch(color='black'    , label=" >= {0}".format(min_library_size))
    lt = mpatches.Patch(color='firebrick', label=" <  {0}".format(min_library_size))
    #plt.legend(handles=[lt, gt], loc='upper center', prop={'size':6}, ncol=2)
    plt.legend(handles=[lt, gt],loc='upper center', bbox_to_anchor=(0.5, -0.12), fancybox=True, ncol=2,  prop={'size':6})

    # Turn of right and top axis
    plt.gca().spines['top'].set_visible(False)     # Turn off top   axis line
    plt.gca().spines['top'].set_color('none')      # Turn off top   axis ticks
    plt.gca().spines['right'].set_visible(False)   # Turn off right axis line
    plt.gca().spines['right'].set_color('none')    # Turn off right axis ticks
    plt.gca().xaxis.set_ticks_position('bottom')
    plt.gca().yaxis.set_ticks_position('left')

    try:
        # Save plot of uniquely mapped reads for every sample
        output_uniquely_mapped_reads_bar_plot = "{0}/{1}_barplot_uniquely_mapped_reads.pdf".format(output_pdf_dir, get_file_info(output_file)[1])
        plt.savefig(output_uniquely_mapped_reads_bar_plot, format="pdf", bbox_inches='tight')
        
        output_uniquely_mapped_reads_bar_plot = get_file_info(output_file)[3] + "_barplot_uniquely_mapped_reads.png"
        plt.savefig(output_uniquely_mapped_reads_bar_plot, format="png", dpi=400, bbox_inches='tight', pad_inches=0.1)

        # clear up the plot
        plt.close('all')
    except Exception as e: 
        print(str(e))
        pass

def get_filtered_data(df, rnaStartCol, rnaEndCol, read_cutoff=5, sample_pc=95, library_size=0):
    ''' filter the data for age and expression '''
    
    # Filter samples with library size
    if library_size > 0:
        df = df[df['uniquely_mapped_reads'] >= library_size]

    # Get the number of colums
    # df.shape gives a tuple with (n_rows, n_columns)
    nrows, ncols = df.shape 

    # Get the miRNA data
    Xrna = df[df.columns[rnaStartCol:ncols-rnaEndCol]] # 2: is to exclude age and gender from this filter
    
    # Filter the reads that are less than say 5
    # find the miRNAs that have 5 or more reads in more than 95% of the data
    try:
        fXrna = Xrna[Xrna<=sp.stats.scoreatpercentile(Xrna, sample_pc, axis=0)].mean(axis=0) >= read_cutoff
        Xrna  = Xrna.loc[:,fXrna]
    except Exception as e: 
        print(str(e))
        pass

    # Add age and gender back
    #Xrna  = pd.concat([df[df.columns[0:rnaStartCol]], Xrna], axis=1) 

    # Get the new dataframe (clinical features, xRNA, remaining clinical features)
    df = pd.concat([df[df.columns[0:rnaStartCol]], Xrna, df[df.columns[ncols-rnaEndCol:ncols]]], axis=1 )

    return df, Xrna

def smallrna_abundance_barcharts(smallrna_type, smallrna_cols, data, group, output_file, n=15):
    ''' Get the pie chart of frequency of mapped reads  '''
    # smallrna_abundance_barcharts(smallrna_type, smallrna_bcol, alzcltdata, 'All_nonCoding_RNA', output_file, 15)

    # Default plot settings
    matplotlib.rcParams['xtick.major.pad']= 15
    matplotlib.rcParams['ytick.major.pad']= 15
    
    for stype, col in zip(smallrna_type, smallrna_cols):
        try:
            print("\t- Plotting {0} ...".format(stype))
            c = Counter()
            for index, row in data.iterrows():
                smrnas = Counter(row[stype])
                c += smrnas
            cd = c.most_common()[:]
            df=pd.DataFrame(cd, dtype='float')
            plt.figure(figsize=(25,15), frameon=False)

            # Get the list of top n(15) labels
            x = df[0]
            y = np.array(df[1])

            # For raw values
            pc   = list(y)
            pcndiff = sum(pc) - sum(pc[:n])
            ypc = pc[:n]
            ypc.append(pcndiff)

            # For percentages
            perc = list(100.*y/y.sum())
            percndiff = sum(perc) - sum(perc[:n])
            yperc = perc[:n]
            yperc.append(percndiff)

            nx  = list(x[:n])
            nx.append("rest (n={0})\nmean = {1:.0f}\nmedian = {2:.0f}".format(len(y)-n, np.mean(pc[n:]), np.median(pc[n:])))

            ind = np.arange(len(ypc))
            width = 0.5

            # Create the bar plot and add labels at the top of each bar
            rects = plt.bar(ind, ypc, width=width, color=col)
            custom_autolabel(rects, plt.gca(), yperc)

            # Aligning rotated xticklabels with their respective xticks: use ha=['right', 'center', 'left']
            plt.xticks(ind + width/2, nx, rotation=45, ha='right', fontsize=20)
            plt.xlim([-0.5,ind.size])

            # Set the fig to be square and set title
            plt.ylabel('Total number of reads log10 (normalized, %)', fontsize=25)
            plt.gca().set_yscale('log', basey=10)
            plt.yticks(fontsize=25)
            plt.title(stype + '\nby frequency of mapped reads in ' + group + ' group')

            # Turn of right and top axis
            plt.gca().spines['top'].set_visible(False)   # Turn off top   axis line
            plt.gca().spines['top'].set_color('none')      # Turn off top   axis ticks
            plt.gca().spines['right'].set_visible(False) # Turn off right axis line
            plt.gca().spines['right'].set_color('none')    # Turn off right axis ticks
            plt.gca().xaxis.set_ticks_position('bottom')
            plt.gca().yaxis.set_ticks_position('left')

            # this is an inset axes over the main axes
            stype_list = list(c.values())
            a = plt.axes([.6, .6, .2, .2])
            vplot = plt.violinplot(stype_list, showmeans=True, showextrema=True, bw_method=0.5)
            for patch in vplot['bodies']:
                patch.set_facecolor("grey")
                patch.set_edgecolor("black")
            mtitle = "Mean = {0:.0f}\nMedian = {1:.0f}".format(np.mean(stype_list), np.median(stype_list))
            plt.title(mtitle, fontsize=16)
            plt.ylabel('Log10(normalized reads)')
            plt.gca().set_yscale('log', basey=10)

            # Save the plot
            output_pdf_plot = "{0}/pdf/{1}_{2}_{3}_barplot_abundance.pdf".format(get_file_info(output_file)[0],get_file_info(output_file)[1], group, stype)
            plt.savefig(output_pdf_plot, format="pdf", bbox_inches='tight')

            output_png_plot = get_file_info(output_file)[3] + group + "_" + stype + "_barplot_abundance"
            plt.savefig(output_png_plot + ".png", format="png", bbox_inches='tight')

            # clear up the plot
            plt.close('all')
        except Exception as e: 
            print(str(e))
            pass

def get_counts(input_dir, allncrna_filename, smallrna_type, m, umrfo):
    ''' Get the counts from counts file for each smallrna class'''
    
    # Get the fileid
    fileid = re.sub('_allncrnaCounts_normalized.txt', '', allncrna_filename)

    if smallrna_type == 'allncrna':
        # Get the raw allncrnaCounts files for library size
        #f = os.popen("find " + input_dir +"/raw_counts/" + smallrna_type + "/ -name *" + str(fileid) + "_allncrnaCounts.txt").read().strip()
        f = "{0}/raw_counts/{1}/{2}_allncrnaCounts.txt".format(input_dir, smallrna_type, fileid)

        # Add the information to m
        library_size = get_smallrna_counts(f, True)
        m[fileid]['uniquely_mapped_reads'] = library_size
        
        # print >> sys.stderr, colored("\t- {0:>9}: ".format(library_size),"green"),           
        # print >> sys.stderr, "{0}".format(fileid)
        print("\t- {0:>9}: {1}".format(library_size, fileid))
        umrfo.write("{0}\t{1}\n".format(fileid, library_size))
    else:
        # Get normalized counts file
        #f = os.popen("find " + input_dir +"/normalized_counts/" + smallrna_type + "/ -name *" + str(fileid) + "_" + smallrna_type + "Counts_normalized.txt").read().strip()
        f = "{0}/normalized_counts/{1}/{2}_{1}Counts_normalized.txt".format(input_dir, smallrna_type, fileid)

        try:
            if os.stat(f).st_size > 0:
                # print >> sys.stderr, colored("\t- Sample_id {0:<5} with counts file:\n\t\t{1}".format(fileid,f), "blue")
                print("\t- Sample_id {0:<5} with counts file:\n\t\t{1}".format(fileid,f))
        
                # Get the smallrna counts dict
                d = get_smallrna_counts(f)

                # Store the information in the super dictionary (m)
                smallrnalist = defaultdict(list) 
                normalized_total_counts = 0
                for smallrna, counts in list(d.items()):
                    normalized_total_counts += counts
                    smallrnalist[smallrna] = counts
                    
                    # add normalized total counts to m
                    #print "normalized_total_counts[{0}]: {1}".format(smallrna_type, normalized_total_counts)
                    m[fileid][smallrna_type + '_total_counts'] = normalized_total_counts

                    # Add the information to m
                    m[fileid][smallrna_type] = smallrnalist
            else:
                print(colored("\t- Sample_id: {0} for {1} \n\t- Counts file is empty...".format(fileid, smallrna_type), "red"))
                m[fileid][smallrna_type] = defaultdict(list)
                m[fileid][smallrna_type + '_total_counts'] = 0
        except Exception as e: 
            print(str(e))
            print(colored("\t- Sample_id: {0} for {1} \n\t- Counts file does not exists...".format(fileid, smallrna_type), "red"))
            m[fileid][smallrna_type] = defaultdict(list)
            m[fileid][smallrna_type + '_total_counts'] = 0

def get_smallrna_counts(input_file, library_size_flag = False):
    ''' Parse the input file and add the information for each miRNA
        If the library size flag is set, return the library size for that file
    '''
    
    # Initialize library size
    library_size = 0

    # Get the information in a dictionary
    d = defaultdict(int)
    with open(input_file,'rU') as fg:
        # head Sample_Blood4mwt_mm_smallrna_sr_Maryam_A_1_mirnaCounts.txt
        # mmu-let-7b-3p	7
        # mmu-let-7b-5p	441
        #....
        # hsa-miR-99b-5p	591
        # no_feature	221439
        # ambiguous	526
        # too_low_aQual	0
        # not_aligned	0
        # alignment_not_unique	1387828
        
        # Get the total read counts
        to_skip_cnt = ['__too_low_Qual', '__not_aligned',  '__alignment_not_unique', 'too_low_Qual', 'not_aligned',  'alignment_not_unique']
        to_skip_all = ['__no_feature', '__ambiguous', 'no_feature', 'ambiguous']
        
        # Loop through rest of the file
        for line in useful_lines(fg):
            featureid = str(line.split("\t")[0])
            counts    = int(line.split("\t")[1])

            # skip unmapped/low quality read counts
            if featureid in to_skip_cnt: continue 
            
            # Get the final counts of all the reads
            library_size += counts
            if featureid in to_skip_all: continue 

            # Add the counts for each featureid
            d[featureid] = counts

    if library_size_flag:
        return library_size
    else:
        return d

def autolabel(rects, ax):
    ''' Attach some text labels '''
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%0.0f'%float(height), ha='center', va='bottom')

def custom_autolabel(rects, ax, custom_list=None):
    # attach some text labels
    for i, rect in enumerate(rects):
        custom_label = height = rect.get_height()
        if custom_list:
            custom_label = "{0} ({1:.0f}%)".format(height, custom_list[i])
        ax.text(rect.get_x() + rect.get_width()/2., 1.02*height, '%s' % custom_label, ha='left', va='bottom', rotation=45, fontsize=20)

# Utility function to report best scores
def report(grid_scores, n_top=3):
    top_scores = sorted(grid_scores, key=itemgetter(1), reverse=True)[:n_top]
    for i, score in enumerate(top_scores):
        print(("Model with rank: {0}".format(i + 1)))
        print(("Mean validation score: {0:.3f} (std: {1:.3f})".format(
              score.mean_validation_score,
              np.std(score.cv_validation_scores))))
        print(("Parameters: {0}".format(score.parameters)))
        print("")

def get_cmap_colors(n, mycmap):
    ''' 
    - From the cmap colors, return a color map of colors
    - http://matplotlib.org/dev/examples/color/colormaps_reference.html
    - All colormaps can be reversed by appending _r. For instance, gray_r is the reverse of gray.
    '''
    import matplotlib.cm as mplcm
    import matplotlib.colors as colors

    NUM_COLORS = n
    cm = plt.get_cmap(mycmap)
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
    colors = sorted([scalarMap.to_rgba(i) for i in range(NUM_COLORS)])
    return colors

def useful_lines(f):
    ''' Filter out useless lines from the blast output file '''
    for line in f:
        line = line.strip()
        if line.startswith('#'):
            continue
        if not line:
            continue
        yield line

def check_options():
    ''' Checks the options to the program '''

    # Create parser object
    parser = argparse.ArgumentParser(add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=textwrap.dedent('''\
        ----------------- SAMPLE USAGE ------------------
        - python scripts/ncRNA_mapping_summary.py -id=input/counts
        - python scripts/ncRNA_mapping_summary.py -id=input/counts -of=output/test.txt
        - python scripts/ncRNA_mapping_summary.py -id=input/counts -of=output/controls/summary_plots/test.txt -sc=test
        - python scripts/ncRNA_mapping_summary.py -id=input/counts -of=output/controls/summary_plots/group_countrols.txt -sc=controls -kc

        -------------------------------------------------
        CONTACT: 
        	Gaurav Jain
        	gaurav.jain@dzne.de
        -------------------------------------------------
        '''))

    # Add arguments 
    parser.add_argument("-id", metavar="--inpdir", help="*Input directory containing the counts file\n - normalized_counts\n - raw_counts: You need to enter the uniquely_mapped_reads_file", dest="input_counts_dir", type=str, required=True)
    parser.add_argument("-of", metavar="--opfile", help="*Output file name"  , dest="output_file" , type=str)
    parser.add_argument("-ls", metavar="--minlib", help="Minimum library size below which plot all the samples red in library size plot.\nDefault = 100000"  , dest="min_library_size" , type=int, default=100000)
    parser.add_argument("-sc", metavar="--subset", help="Comma separated file name regexes to plot only subset of these files\n- Examples:\n\t1) _A_\n\t2) NDC\n\t3) '_B_|_C_'\n\t4) Exo.*V1", dest="subset_files", type=str)
    parser.add_argument("-sf", metavar="--sbfile", help="Subset filenames list. For exmaple file containing AD samples:\n\t119_hs_smallrna_sr_Pooja_C_1_1\n\t123_hs_smallrna_sr_Pooja_C_2_1\n\t144_hs_smallrna_sr_Pooja_C_4_1"  , dest="subset_sample_file" , type=str)
    parser.add_argument('-mx', "--merged_excel"  , help="if set, also make the merged excel of all smallrnas summary\nNOTE: Takes bit long for 20 or more samples...", action='store_true', dest="merged_excel", default=False)
    parser.add_argument('-kc', "--skip_conf"     , help="if set, skip conformation from user for the subset datasets", action='store_true', dest="skip_subset_conf", default=False)

    # Print the help message only if no arguments are supplied
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    # Save the STDOUT output in a log file
    if parser.parse_args().output_file:
        logdir = "{0}/logs".format(get_file_info(parser.parse_args().output_file)[0])
        create_dir(logdir)
        logfile = "{0}/{1}.log".format(logdir, get_file_info(parser.parse_args().output_file)[1])
    else:
        #logdir  = "{0}/logs".format(os.getcwd())
        logdir  = "{0}/summary_plots/all/logs".format(parser.parse_args().input_counts_dir)
        create_dir(logdir)
        logfile = "{0}/{1}.log".format(logdir,get_file_info(sys.argv[0])[1])

    logf = open(logfile, 'w')
    sys.stdout = Log(logf, sys.stdout)

    # Parse command line with parse_args and store it in an object
    args = parser.parse_args()
    print_initial_arguments(parser)
    return args

# main function
if __name__=="__main__":
      main()

