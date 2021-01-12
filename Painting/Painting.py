# Not finished version of painting

# TODO
# 1) Warp in functions
# 2) Usage from command line
# 3) Current programme needs a lot of memory. Try to use libraries for reading the data (pandas/numpy)
# 4) Try seaborn, plotly and bokeh for visualisation
# 5) Use another format of data (fasta), not only fastq
# 6) Add data name to names of paintings (??)

from collections import Counter
import statistics as stat
import numpy as np
import os

import matplotlib.pyplot as plt
import seaborn as sns
import plotly
import bokeh

# 1. Per base sequence quality
# 2. Per sequence quality scores
# 3. Per base sequence content
# 4. Per sequence GC content
# 5. Per base N content
# 6. Sequence Length Distribution

read_quality_per_base_huge = [[] for _ in range(1000)]  # 1. Per base sequence quality
read_quality_per_base_mean = []  # 1. Per base sequence quality

read_quality_per_seq = []  # 2. Per sequence quality scores

sequence_content_huge = [[] for _ in range(1000)]  # 3. Per base sequence content and 5. Per base N content
A_per_base = []  # 3. Per base sequence content
T_per_base = []  # 3. Per base sequence content
G_per_base = []  # 3. Per base sequence content
C_per_base = []  # 3. Per base sequence content
GC_content = []  # 4. Per sequence GC content
N_per_base = []  # 5. Per base N content
read_len = []  # 6. Sequence Length Distribution

# with open("Raddei.fastq", "r") as fastq_data:  # 5.6 Gb of data. Not working !!!
with open("Test_reads_full_by_AnyaChurkina.fastq", "r") as fastq_data:  # Full
# with open("Test_reads_cut_by_AnyaChurkina.fastq", "r") as fastq_data:  # Cut (Len of reads < 50 p.b.)
    print("Start reading")
    short_data_counter = 0  # For testing

    line_number = 0
    for line in fastq_data:
        short_data_counter += 1  # For testing
        line_number += 1

        if line_number == 2:  # String with nucleotides

            # Count GC (4. Per sequence GC content)
            GC = 100 * ((line.rstrip().count("G")) + (line.rstrip().count("C"))) / len(line.rstrip())
            GC_content.append(int(GC))

            # For T, G, C, A and N content (3. Per base sequence content and 5. Per base N content)
            for i in range(len(line.rstrip())):
                sequence_content_huge[i].append(line.rstrip()[i])

            # For len distribution (6. Sequence Length Distribution)
            read_len.append(len(line.rstrip()))

        elif line_number == 4:  # String with quality

            # Get quality (Phred Score)
            phred_value = [ord(ch) - 33 for ch in line.rstrip()]

            # For "1. Per base sequence quality"
            for i in range(len(phred_value)):
                read_quality_per_base_huge[i].append(phred_value[i])

            # For "2. Per sequence quality scores"
            read_quality_per_seq.append(int(stat.mean(phred_value)))

            line_number = 0

        # For testing _______
        # if short_data_counter > 64:
        #     break

print("Reading data is done")


# ______________ Preprocessing ______________
print("Start processing the data")
# 1. Per base sequence quality
read_quality_per_base = read_quality_per_base_huge[:read_quality_per_base_huge.index([])]

if len(read_quality_per_base) > 50:
    read_quality_per_base_cut = read_quality_per_base[:9]
    xlab_per_base_seq_qual = [str(i) for i in range(1, 10)]

    for i in range(15, len(read_quality_per_base), 5):
        sub_list = []
        for sub in read_quality_per_base[i - 5:i - 1]:
            for item in sub:
                sub_list.append(item)

        read_quality_per_base_cut.append(sub_list)
        xlab_per_base_seq_qual.append(f"{i - 5}-{i - 1}")

else:
    read_quality_per_base_cut = read_quality_per_base[:]
    xlab_per_base_seq_qual = [str(i) for i in range(1, len(read_quality_per_base) + 1)]

for val_list in read_quality_per_base_cut:
    read_quality_per_base_mean.append(int(stat.mean(val_list)))

x_qpb = [i for i in range(len(xlab_per_base_seq_qual))]

# xlab_per_base_seq_qual - name of each x value
# x_qpb - x value
# read_quality_per_base_cut - y values for boxplots (for each x value)
# read_quality_per_base_mean - y values for mean (for each x value)


# 2. Per sequence quality scores
sorted_quality_data = sorted(read_quality_per_seq)
count_quality = Counter(sorted_quality_data)

x_per_seq_qual_scores = list(count_quality.keys())  # x value
y_per_seq_qual_scores = list(count_quality.values())  # y value


# 3. Per base sequence content and 5. Per base N content
sequence_content = sequence_content_huge[:sequence_content_huge.index([])]

if len(sequence_content) > 50:
    sequence_content_cut = sequence_content[:19]
    xlab_seq_cont = [str(i) for i in range(1, 20)]

    for i in range(25, len(sequence_content), 5):
        sub_list = []
        for sub in sequence_content[i - 5:i - 1]:
            for item in sub:
                sub_list.append(item)

        sequence_content_cut.append(sub_list)
        xlab_seq_cont.append(f"{i - 5}-{i - 1}")

else:
    sequence_content_cut = sequence_content[:]
    xlab_seq_cont = [str(i) for i in range(1, len(sequence_content) + 1)]


for content_per_base in sequence_content_cut:
    A_val = 100 * (content_per_base.count("A") / len(content_per_base))
    T_val = 100 * (content_per_base.count("T") / len(content_per_base))
    G_val = 100 * (content_per_base.count("G") / len(content_per_base))
    C_val = 100 * (content_per_base.count("C") / len(content_per_base))
    N_val = 100 * (content_per_base.count("N") / len(content_per_base))

    A_per_base.append(A_val)  # y value
    T_per_base.append(T_val)  # y value
    G_per_base.append(G_val)  # y value
    C_per_base.append(C_val)  # y value
    N_per_base.append(N_val)  # y value

    # xlab_seq_cont - x value


# 4. Per sequence GC content
GC_content = sorted(GC_content)
count_GC_content = Counter(GC_content)

x_GC = list(count_GC_content.keys())  # x value
y_GC = list(count_GC_content.values())  # y value


# 6. Sequence Length Distribution
sorted_len_data = sorted(read_len)
# count = {i: len_data.count(i) for i in len_data_sort}
count_len = Counter(sorted_len_data)

x_len = list(count_len.keys())  # x value
y_len = list(count_len.values())  # y value

print("Preprocessing is done")


# ______________ Directory for paintings ______________

path = "./Matplotlib_paintings"

if not os.path.exists(path):
    os.mkdir(path)
    print(f"Creating {os.getcwd() + path[1:]}")

# ______________ PLOTS ______________

# 1. Per base sequence quality
fig, ax = plt.subplots(figsize=(20, 10))

bp = ax.boxplot(read_quality_per_base_cut, positions=x_qpb, patch_artist=True,
                widths=0.7, showfliers=False, labels=xlab_per_base_seq_qual)

ax.plot(x_qpb, read_quality_per_base_mean, "b-")

fig.autofmt_xdate()

for patch in bp['boxes']:
    patch.set(facecolor="yellow")
    patch.set(alpha=0.7)

start, end = ax.get_ylim()
ax.yaxis.set_ticks(np.arange(int(0), int(end), 2))

ax.fill_between(x_qpb, 0, 20, alpha=0.5, color="lightcoral")
ax.fill_between(x_qpb, 20, 28, color="lightyellow")
ax.fill_between(x_qpb, 28, 40, alpha=0.5, color="lightgreen")


plt.xlabel("Position in read (bp)")
plt.ylabel("Quality (Phred Score)")
plt.title("Quality scores across all bases")
plt.grid(axis="y")

plt.savefig(path + "/1_Per_base_sequence_quality.jpg", format="jpg", bb_inches="tight")
plt.close()


# 2. Per sequence quality scores
plt.figure(figsize=(20, 10))
plt.plot(x_per_seq_qual_scores, y_per_seq_qual_scores, color="red")
plt.xlabel("Mean Sequence Quality (Phred Score)")
plt.ylabel("Number of Sequence")
plt.title("Quality score distribution over all sequences")
plt.grid(axis="y")

plt.savefig(path + "/2_Per_sequence_quality_scores.jpg", format="jpg", bb_inches="tight")
plt.close()


# 3. Per base sequence content
fig, ax = plt.subplots(figsize=(20, 10))

ax.plot(xlab_seq_cont, T_per_base, color="red")
ax.plot(xlab_seq_cont, C_per_base, color="blue")
ax.plot(xlab_seq_cont, A_per_base, color="green")
ax.plot(xlab_seq_cont, G_per_base, color="black")

plt.xlabel("Position in read (bp)")
plt.gca().set_ylim([0, 100])
plt.title("Sequence content across all bases")
plt.grid(axis="y")
fig.autofmt_xdate()
plt.legend(("%T", "%C", "%A", "%G"), loc="upper right")

plt.savefig(path + "/3_Per_base_sequence_content.jpg", format="jpg", bb_inches="tight")
plt.close()


# 4. Per sequence GC content
plt.figure(figsize=(20, 10))
plt.plot(x_GC, y_GC, color="red")
plt.xlabel("Mean GC content (%)")
plt.ylabel("Number of Sequence")
plt.title("GC content distribution over all sequences")
plt.grid(axis="y")

plt.savefig(path + "/4_Per_sequence_GC_content.jpg", format="jpg", bb_inches="tight")
plt.close()


# 5. Per base N content
fig, ax = plt.subplots(figsize=(20, 10))
ax.plot(xlab_seq_cont, N_per_base, color="red")

plt.xlabel("Position in read (bp)")
plt.gca().set_ylim([0, 100])
plt.title("N content across all bases")
plt.grid(axis="y")
fig.autofmt_xdate()

plt.savefig(path + "/5_Per_base_N_content.jpg", format="jpg", bb_inches="tight")
plt.close()


# 6. Sequence Length Distribution
plt.figure(figsize=(20, 10))
plt.plot(x_len, y_len, color="red")

plt.xlabel("Sequence Length (bp)")
plt.ylabel("Number of Sequence")
plt.title("Sequence Length Distribution")
plt.grid(axis="y")

plt.savefig(path + "/6_Sequence_Length_Distribution.jpg", format="jpg", bb_inches="tight")
plt.close()

print()
print("Painting with Matplotlib is done")
print(f"Paintings were placed in {os.getcwd() + path[1:]}")
