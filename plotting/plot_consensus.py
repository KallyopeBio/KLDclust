import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
from intervaltree import Interval, IntervalTree

## Things to try
## Given # of interval, what is the min % excluded
## Given max length of an interval, what is the min % excluded
## Given the % coverage desired, what is the max # intervals
## Given the % coverage desired, what is the min length of an interval
## Weight each variant by importance, minimize weight excluded

## SET INPUT DATA
## Only parameter impacting the consensus analysis in this code
buffer_len = 5000

## Provide links to data
MAF="05"
root_dir="/home/sarah/Data/1000genomes_phase3/Fisher_Penalty/MAF"+MAF+"/chr21/"
file_name="clustering/FINAL_optBreaks.txt"
AFR_file=root_dir+"AFR/"+file_name
AMR_file=root_dir+"AMR/"+file_name
EAS_file=root_dir+"EAS/"+file_name
EUR_file=root_dir+"EUR/"+file_name
SAS_file=root_dir+"SAS/"+file_name

## Set information for chart visual and title
lower_bound=22000000
upper_bound=42000000
y_bound=2000000
file_list=[AFR_file, SAS_file, EAS_file, EUR_file]
pop_list=['AFR', 'SAS', 'EAS', 'EUR']
color_list=['r', 'g', 'c', 'b']


## DEFINE FUNCTIONS
def get_starts(input_file):
    start_list = []
    with open(file, "r") as fin:       
        for line in fin:
            line.strip("\n")
            endpoints=line.split(" ")
            start = int(endpoints[0])
            start_list.append(start)
    return start_list

def get_stops(input_file):
    stop_list = []
    with open(file, "r") as fin:       
        for line in fin:
            line.strip("\n")
            endpoints=line.split(" ")
            stop = int(endpoints[1])
            stop_list.append(stop)
    return stop_list

def get_extremes(start_points, stop_points):
    chrom_start=float('inf')
    chrom_stop=0

    for key, value in start_data.iteritems():
        for i in range(0, len(value)):
            if(value[i]<chrom_start):
                chrom_start=value[i]
            if(stop_data[key][i]>chrom_stop):
                chrom_stop=stop_data[key][i]
    return (chrom_start, chrom_stop)
           
def get_overlap(interval_1, interval_2):
    start_1=interval_1[0]
    start_2=interval_2[0]
    stop_1=interval_1[1]
    stop_2=interval_2[1]

    if start_1>start_2:
        overlap_start=start_1
    else:
        overlap_start=start_2
    if stop_1<stop_2:
        overlap_stop=stop_1
    else:
        overlap_stop=stop_2
    return (overlap_start, overlap_stop)

def plot_intervals(start_points, stop_points, color, label, linestyle='solid'):
    assert(len(start_points)==len(stop_points))
    x=[]
    y=[]
    for i in range(0, len(start_points)):
        start = start_points[i]
        stop = stop_points[i]
        height = (stop-start)/2.0
        midpoint = start + height
        x.extend((start, midpoint, stop))
        y.extend((0, height, 0))

    print(label, len(x)/3)
    plt.plot(x, y, color=color, label=label, linewidth=.5, linestyle=linestyle)

    return
   
def interval_tree(start_data, stop_data, buffer_len):
    starts=[]
    stops=[]
    t = IntervalTree()

    ## Shrink each interval by the buffer size
    for key, value in start_data.iteritems():
        for i in range(0, len(value)):
            shrunk_start=value[i]+buffer_len/2.0
            shrunk_stop=stop_data[key][i]+1-buffer_len/2.0
            if shrunk_start < shrunk_stop:
                t[shrunk_start:shrunk_stop]=(shrunk_start, shrunk_stop)

    ## Add chromosome endpoints without buffer 
    chrom_start, chrom_stop = get_extremes(start_data, stop_data)
    if chrom_start<t.begin()+1:
        t[chrom_start:t.begin()+1]=(chrom_start, t.begin()+1)
    if t.end()-1<chrom_stop:
        t[t.end()-1:chrom_stop]=(t.end()-1, chrom_stop)

    ## Merge intervals that overlap in tree to get consensus
    t.merge_overlaps()

    ## Check that original intervals only overlap with one consensus interval
    for key, value in start_data.iteritems():
        for i in range(0, len(value)):
            start=value[i]
            stop=stop_data[key][i]+1
            if len(t[start:stop])>1:
                ## If they overlap with more than one
                ## Remove part of consensus interval
                ## This will never be more than the buffer size/2
                assert(len(t[start:stop])==2)
                remove_start=0
                remove_stop=0
                min_length=float('inf')
                for interval in t[start:stop]:
                    overlap_start, overlap_stop=get_overlap((start, stop), (interval[0], interval[1]))
                    if (overlap_stop-overlap_start)<min_length:
                        min_length=overlap_stop-overlap_start
                        remove_start=overlap_start
                        remove_stop=overlap_stop
                print(min_length)
                t.chop(remove_start, remove_stop)
                assert(min_length<=buffer_len/2.0)
                assert(len(t[start:stop])<2)
    
    ## Get consensus start and stop points
    chrom_len=chrom_stop-chrom_start
    covered=0.0
    for interval in sorted(t):
        starts.append(interval[0])
        stops.append(interval[1])
        covered=covered+(interval[1]-interval[0])

    print("The percentage of the chromosome covered is: %s" % '{0:.2f}'.format((covered/chrom_len)*100.0))
        
    return (starts, stops)
    
## MAIN FUNCTION
i=0
start_data={}
stop_data={}

for file in file_list:

    ## Add data (without midpoints) to dictionary
    start_data[file]= get_starts(file)
    stop_data[file]=get_stops(file)
    
    ## Now we add data to plot
    plot_intervals(get_starts(file), get_stops(file), color=color_list[i], label=pop_list[i])
    i = i+1

## Make interval tree
consensus_starts, consensus_stops = interval_tree(start_data, stop_data, buffer_len)
plot_intervals(consensus_starts, consensus_stops, 'y', 'Consensus', 'dashed')
    
## Now that data has been processed we format and plot
plt.gca().xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
plt.xlabel('Chromosome (BP)')
plt.xlim((lower_bound, upper_bound))
plt.ylim((0,y_bound))
plt.title('LD Blocks by Population with MAF%s' %MAF)
plt.legend()
plt.rcParams["figure.figsize"] = (10,5)
plt.show()
plt.gcf().clear()
