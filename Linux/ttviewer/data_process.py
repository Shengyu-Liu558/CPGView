import re
import os
import csv
import sys
from Bio import SeqIO

# Part one: import data
arch = sys.argv[1]
file1 = sys.argv[2]
file2 = sys.argv[3]


def get_path(path1):
    r = os.path.abspath(path1)
    return r


def overlap(new_list22):
    new_list = []
    new_list1 = []
    for aq in new_list22:
        qs1 = [aq[0], aq[1], aq[2], aq[3], aq[4], aq[5], aq[6], aq[7], aq[6], aq[7], aq[-1]]
        new_list.append(qs1)

    m = {}
    for i in new_list:

        key = i[0]
        if key in m:
            m[key].append(i[0:])

        else:
            m[key] = [i[0:]]

    for j in m:
        total_distancem = (m[j][0][2] - m[j][0][1]) / 8
        total_distancen = m[j][0][2] - m[j][0][1]
        he_list = []
        new_he_list = []

        for s in m[j]:
            he_list.append(s[-5])
            he_list.append(s[-4])

        for sj in range(0, len(he_list) - 1):

            mj = sj + 1
            if mj != len(he_list) - 1:
                if he_list[mj] - he_list[sj] < total_distancem:

                    res = he_list[mj] - he_list[sj]
                    he_list[mj] = he_list[mj] - res + total_distancem

                else:
                    continue
            else:
                if he_list[mj] >= he_list[sj]:
                    res1 = he_list[mj] - he_list[sj]
                    he_list[mj] = he_list[mj] - res1 + total_distancem
                else:
                    res2 = he_list[sj] - he_list[mj]
                    he_list[mj] = he_list[mj] + res2 + total_distancem

        new_distance = he_list[-1] - he_list[0]

        for yu in range(1, len(he_list)):
            he_list[yu] = int((((he_list[yu] - he_list[0]) * total_distancen) / new_distance) + m[j][0][1])

        if abs(he_list[0] - m[j][0][1]) < 2:
            he_list[0] = m[j][0][1]
        if abs(he_list[-1] - m[j][0][2]) < 2:
            he_list[-1] = m[j][0][2]

        for ms1 in range(0, len(he_list), 2):
            l2 = he_list[ms1: ms1 + 2]
            new_he_list.append(l2)

        for ye1, ye2 in zip(m[j], new_he_list):
            ye1[-5] = ye2[0]
            ye1[-4] = ye2[1]

        for js in m[j]:
            new_list1.append(js)

    return new_list1


def add_intron(new_list23):
    test_value1 = new_list23[0]
    test_value2 = new_list23[1]

    if test_value1[0] == test_value2[0] and test_value1[1] == test_value2[1] and test_value1[2] == test_value2[2]:
        intron = [test_value1[0], test_value1[1], test_value1[2], test_value1[3], test_value1[4], "Intron",
                  test_value1[7], test_value2[6], test_value1[9], test_value2[8], test_value1[-1]]

    else:
        if test_value1[1] == test_value1[-2]:
            intron = [test_value1[0], test_value1[1], test_value1[2], test_value1[3], test_value1[4], "Intron",
                      test_value1[7], test_value1[2], test_value1[9], test_value1[2], test_value1[-1]]
        elif test_value1[2] == test_value1[-1]:
            intron = [test_value1[0], test_value1[1], test_value1[2], test_value1[3], test_value1[4], "Intron",
                      test_value1[1], test_value1[6], test_value1[1], test_value1[8], test_value1[-1]]

    new_list23.insert(0, intron)

    return new_list23


def check_rps12(new):
    M = []
    for i in new:
        if i[2] < i[1]:

            M.append(1)
        else:
            continue
    return M


def mkdir(path):
    folder = os.path.exists(path)
    if not folder:
        os.makedirs(path)
        print("")
        print("---  New folder (", str(path).strip(), ") has been generated  ---")
        return 1
    else:
        print("")
        print("---  New folder (", str(path).strip(), ") already exists here. Please delete it first!  ---")
        return -1


error0 = ""
try:
    record = SeqIO.parse(arch, "genbank")
    rec = next(record)
    for rec in SeqIO.parse(arch, 'gb'):
        print(rec.annotations)

    new_list1 = []
    new_list2 = []
    new_list3 = []
    new_list4 = []
    new_list5 = []
    new_list6 = []

    # Part two: select gene and subgene
    y_nums = 0
    for f in rec.features:
        direction = 0
        strand = ""

        if f.type != 'source' and f.type.lower() == 'cds' and 'gene' in f.qualifiers:
            if "join" in str(f.location) and "('gene'," in str(f.qualifiers):

                y_nums += 1
                location1 = str(f.location).replace("<", "").replace(">", "").replace("join{", "").replace("}",
                                                                                                           "").replace(
                    "](-)", "").replace("](+)", "").replace("[", "").replace(" ", "").split(
                    ",")  # Remove the position of the direction
                sign = list(set(re.findall(r'[(](.*?)[)]', str(f.location))))[0]

                if sign == "-":
                    direction = -1
                    strand = "reverse"
                elif sign == "+":
                    direction = 1
                    strand = "forward"

                if "(+)" in str(f.location):

                    if f.qualifiers["gene"][0].lower() != 'rps12':
                        location1 = sorted(location1, key=lambda i: int(i.split(":")[0]), reverse=False)
                    else:
                        location1 = location1
                    start1 = int(location1[0].split(":")[0]) + 1
                    end1 = int(location1[-1].split(":")[1])
                    line = [str(y_nums), int(start1), int(end1), strand, direction, f.qualifiers["gene"][0]]
                    new_list1.append(line)

                    for jj in location1:
                        substart1 = int(jj.split(":")[0]) + 1
                        subend1 = jj.split(":")[1]
                        line3 = [str(y_nums), int(start1), int(end1), strand, direction, "Exon", int(substart1),
                                 int(subend1), f.qualifiers["gene"][0]]
                        new_list2.append(line3)

                elif "(-)" in str(f.location):

                    if f.qualifiers["gene"][0].lower() != 'rps12':
                        location1 = sorted(location1, key=lambda i: int(i.split(":")[0]), reverse=True)
                    else:
                        location1 = location1

                    start2 = int(location1[-1].split(":")[0]) + 1
                    end2 = int(location1[0].split(":")[1])
                    line1 = [str(y_nums), int(start2), int(end2), strand, direction, f.qualifiers["gene"][0]]
                    new_list1.append(line1)

                    for jj in location1:
                        substart2 = int(jj.split(":")[0]) + 1
                        subend2 = jj.split(":")[1]

                        line2 = [str(y_nums), int(start2), int(end2), strand, direction, "Exon", int(substart2),
                                 int(subend2), f.qualifiers["gene"][0]]
                        new_list2.append(line2)

                if f.qualifiers["gene"][0].lower() == 'rps12':

                    location2 = str(f.location).replace("<", "").replace(">", "").replace("join{", "").replace("}",
                                                                                                               "").replace(
                        "]", "").replace("]",
                                         "").replace(
                        "[", "").replace(" ", "").split(",")  # Remove the position of the direction
                    sign2 = list(set(re.findall(r'[(](.*?)[)]', str(f.location))))
                    unms = 1

                    new_location = []
                    for mie in location2:
                        mie1 = mie + str(unms)
                        new_location.append(mie1)
                        unms += 1

                    str_new_location = " ".join(new_location)

                    if "+" in str_new_location and "-" in str_new_location:

                        new_location1 = str_new_location.split(" ")
                        for mxs in new_location1:
                            if "+" in mxs:
                                direction1 = 1
                                strand1 = "forward"
                                mx1 = mxs.split("(")

                                label_start1 = int(mx1[0].split(":")[0]) + 1
                                label_end1 = int(mx1[0].split(":")[1])
                                num1 = mx1[1].split(")")[1]

                                line3 = ["2 rps12", "exon" + str(num1), strand1, direction1, int(label_start1),
                                         int(label_end1)]
                                new_list3.append(line3)

                                new_line3 = ["1 Transcript 1", "exon" + str(num1), "Transcript 1", direction1,
                                             int(label_start1), int(label_end1), "exon" + str(num1)]

                                if new_line3[1] == "exon2" or new_line3[1] == "exon3":
                                    new_line3[1] = new_line3[1] + " (IRa)"
                                new_list4.append(new_line3)

                            elif "-" in mxs:
                                direction2 = -1
                                strand2 = "reverse"
                                mx2 = mxs.split("(")
                                label_start2 = int(mx2[0].split(":")[0]) + 1
                                label_end2 = int(mx2[0].split(":")[1])
                                num2 = mx2[1].split(")")[1]
                                line4 = ["2 rps12", "exon" + str(num2), strand2, direction2, int(label_start2),
                                         int(label_end2)]
                                new_list3.append(line4)

                                new_line4 = ["1 Transcript 1", "exon" + str(num2), "Transcript 1", 1, int(label_start2),
                                             int(label_end2), "exon" + str(num2)]
                                new_list4.append(new_line4)

                    elif "-" in str_new_location and "+" not in str_new_location:

                        new_location2 = str_new_location.split(" ")
                        for mxs1 in new_location2:
                            direction3 = -1
                            strand3 = "reverse"
                            mx3 = mxs1.split("(")
                            label_start3 = int(mx3[0].split(":")[0]) + 1
                            label_end3 = int(mx3[0].split(":")[1])
                            num3 = mx3[1].split(")")[1]
                            line5 = ["2 rps12", "exon" + str(num3), strand3, direction3, int(label_start3),
                                     int(label_end3)]
                            new_list3.append(line5)

                            new_line5 = ["3 Transcript 2", "exon" + str(num3), "Transcript 2", 1, int(label_start3),
                                         int(label_end3), "exon" + str(num3)]
                            if new_line5[1] == "exon2" or new_line5[1] == "exon3":
                                new_line5[1] = new_line5[1] + " (IRb)"
                            new_list5.append(new_line5)

                    elif "+" in str_new_location and "-" not in str_new_location:
                        new_location3 = str_new_location.split(" ")
                        for mxs2 in new_location3:
                            direction4 = 1
                            strand4 = "forward"
                            mx4 = mxs2.split("(")
                            label_start4 = int(mx4[0].split(":")[0]) + 1
                            label_end4 = int(mx4[0].split(":")[1])
                            num4 = mx4[1].split(")")[1]
                            line6 = ["2 rps12", "exon" + str(num4), strand4, direction4, int(label_start4),
                                     int(label_end4)]
                            new_list3.append(line6)
                            new_line6 = ["3 Transcript 2", "exon" + str(num4), "Transcript 2", 1, int(label_start4),
                                         int(label_end4), "exon" + str(num4)]
                            if new_line6[1] == "exon2" or new_line6[1] == "exon3":
                                new_line6[1] = new_line6[1] + " (IRa)"
                            new_list6.append(new_line6)

    new_list5 = (new_list6 if new_list6 != [] else new_list5)

    error3x = []
    if new_list1 == []:
        error1 = "<h4>Error!</h4><p>There are no gene information or there is gene information, but there is no CDS information in the GenBank file.</p> "

    else:
        error1 = ""
        for check_error1 in new_list1:

            if "rps12" == check_error1[-1]:
                error3x.append(0)
                continue
            else:
                error3x.append(1)

    if 0 not in error3x:
        error3 = "<h4>Error!</h4><p>No rps12 gene is found in GenBank file. Alternatively, the rps12 gene is found, but no introns are found for the rps12 gene.</p> "
    else:
        error3 = ""

    m_numb1 = 0
    for woe in new_list1:
        m_numb1 = m_numb1 + (woe[5].lower().count("rps12"))

    check_gene = check_rps12(new_list1)

    for index_gene0 in new_list1:
        if index_gene0[-1].lower() == 'rps12':
            new_list1.remove(index_gene0)
    for index_gene1 in new_list2:
        if index_gene1[-1].lower() == 'rps12':
            new_list2.remove(index_gene1)

    # Part three: delete  genes and errors in genes
    new_list1 = sorted(new_list1, key=lambda x: x[1])
    new_list11 = new_list1.copy()
    num_list = []

    for i in range(0, len(new_list1) - 1):
        for ii in range(i + 1, len(new_list1)):
            if new_list1[i][2] > new_list1[ii][1]:
                if new_list1[i] in new_list11:
                    new_list11.remove(new_list1[i])
                else:
                    break
            else:
                continue

    for ww1 in range(0, len(new_list11) - 1):
        if ww1 != len(new_list11) - 1:
            ww2 = ww1 + 1
            if new_list11[ww2][1] - new_list11[ww1][2] < 300:
                rwqs = 300 - (new_list11[ww2][1] - new_list11[ww1][2])
                for ww3 in range(ww2, len(new_list11)):
                    new_list11[ww3][1] = new_list11[ww3][1] + rwqs
                    new_list11[ww3][2] = new_list11[ww3][2] + rwqs
        else:
            break

    for rr in new_list11:
        num_list.append(rr[0])

    if new_list11 == []:
        error2 = "<h4>Error!</h4><p>The genes contained in the GenBank file are not arranged in the correct order, i.e., from small to large.</p> "
    else:
        error2 = ""

    new_list2 = sorted(new_list2, key=lambda x: x[6])

    new_list22 = new_list2.copy()

    for m in new_list2:
        if m[0] not in num_list:
            new_list22.remove(m)

        else:
            continue

    # Part four: Dealing with the problem of gene exon coordinate overlap

    new_list23 = overlap(new_list22)

    for jw1 in new_list23:
        for jw2 in new_list11:
            if jw1[0] == jw2[0] and jw1[-1] == jw2[-1]:

                chax1 = jw1[2] - jw1[1]

                chax2 = jw1[6] - jw1[1]
                chax3 = jw1[7] - jw1[1]

                bili1 = chax2 / chax1
                bili2 = chax3 / chax1

                jw1[6] = int((bili1 * chax1) + jw2[1])
                jw1[7] = int((bili2 * chax1) + jw2[1])

                jw1[1] = jw2[1]
                jw1[2] = jw2[2]

                if abs(jw1[6] - jw1[1]) < 2:
                    jw1[6] = jw1[1]
                if abs(jw1[7] - jw1[2]) < 2:
                    jw1[7] = jw1[2]
    if new_list23 != []:
        new_list23 = add_intron(new_list23)
    else:
        new_list23 = new_list23

    # The sixth part: processing the data of rps12

    new_raw_list0 = []
    new_raw_list1 = []
    new_raw_list2 = []
    new_list33 = new_list3.copy()
    list1 = list(set(list(map(tuple, new_list33))))
    list2 = list(map(list, list1))
    list2 = sorted(list2, key=lambda x: x[-2])
    list3 = []
    list4 = []
    for lw in list2:
        if lw[3] == 1:
            list3.append(lw)
        elif lw[3] == -1:
            list4.append(lw)

    list3 = sorted(list3, key=lambda x: x[-2])  # +
    list4 = sorted(list4, key=lambda x: x[-2])  # -
    check_len = len(list3) + len(list4)
    if list3 != [] and list4 != []:

        local1 = []
        local2 = []
        local3 = []
        local4 = []
        new_raw0 = []
        new_trans0 = []
        if check_len == 5:
            local1 = [[1260, 1410], [1490, 1640]]
            local2 = [[600, 750], [820, 970], [1040, 1190]]
            local3 = [[100, 135], [255, 290], [310, 345]]
            local4 = [[2100, 2135], [2210, 2245], [2310, 2345]]
        elif check_len == 3:
            local1 = [[1250, 1500], [1490, 1640]]
            local2 = [[500, 750], [870, 1120]]
            local3 = [[100, 150], [250, 300]]
            local4 = [[2100, 2150], [2250, 2300]]

        elif check_len == 2 and (new_list4 == [] or new_list5 == []):
            if new_list3[-1][-2] > new_list3[0][-2]:
                new_raw0 = [["2 rps12", "exon2 (IRLC)", "", "", "Genome  (+)", 1, "", "", "exon2"],
                            ["2 rps12", "exon1 (IRLC)", 500, 750, "(-)", -1, 69007, 69120, "exon1"],
                            ["2 rps12", "exon2 (IRLC)", 870, 1120, "(-)", -1, 96927, 97169, "exon2"]]
                new_trans0 = [["3 Transcript 2", "exon1 (IRLC)", 2100, 2350, "Transcript", 1, 69007, 69120, "exon1"],
                              ["3 Transcript 2", "exon2 (IRLC)", 2470, 2720, "Transcript", 1, 96927, 97169, "exon2"]]

            elif new_list3[-1][-2] < new_list3[0][-2]:
                new_raw0 = [["2 rps12", "exon2 (IRLC)", "", "", "Genome  (+)", 1, "", "", "exon2"],
                            ["2 rps12", "exon1 (IRLC)", "500,750,(-)", 1, 69007, 69120, "exon1"],
                            ["2 rps12", "exon2 (IRLC)", 870, 1120, "(-)", 1, 96927, 97169, "exon2"]]
                new_trans0 = [["3 Transcript 2", "exon1 (IRLC)", 2100, 2350, "Transcript", 1, 69007, 69120, "exon1"],
                              ["3 Transcript 2", "exon2 (IRLC)", 2470, 2720, "Transcript", 1, 96927, 97169, "exon2"]]
            else:
                new_raw0 = []
                new_trans = []

        if local1 == [] and local2 == [] and local3 == [] and local4 == []:
            new_raw_list0 = []
            new_raw_list1 = []
            new_raw_list2 = []

        else:
            for xq, yq in zip(list3, local1):
                xq.insert(2, yq[0])
                xq.insert(3, yq[1])
                xq.insert(8, xq[1])

                xq[4] = "Genome  (+)"
                if xq[1] == "exon2":
                    xq[1] = "exon2" + " (IRa)"
                elif xq[1] == "exon3":
                    xq[1] = "exon3" + " (IRa)"
                new_raw_list1.append(xq)

            for xq1, yq1 in zip(list4, local2):
                xq1.insert(2, yq1[0])
                xq1.insert(3, yq1[1])
                xq1.insert(8, xq1[1])
                xq1[4] = "(-)"
                if xq1[1] == "exon2":
                    xq1[1] = "exon2" + " (IRb)"
                elif xq1[1] == "exon3":
                    xq1[1] = "exon3" + " (IRb)"
                new_raw_list1.append(xq1)

            for xq2, yq2 in zip(new_list4, local3):
                xq2.insert(2, yq2[0])
                xq2.insert(3, yq2[1])
                new_raw_list0.append(xq2)

            for xq3, yq3 in zip(new_list5, local4):
                xq3.insert(2, yq3[0])
                xq3.insert(3, yq3[1])
                new_raw_list2.append(xq3)

    elif (list3 == [] and list4 != []) or (list3 != [] and list4 == []):
        loca21 = []
        if check_len == 3:
            loca21 = [[500, 750], [820, 1120], [1250, 1500]]

            if loca21 == []:
                new_raw_list1 = []
            else:
                if list3 != []:
                    for op, sp in zip(list3, loca21):
                        op.insert(2, sp[0])
                        op.insert(3, sp[1])
                        op.insert(8, op[1])
                        op[4] = "Genome  (+)"

                        new_raw_list1.append(op)
                if list4 != []:
                    for op1, sp1 in zip(list4, loca21):
                        op1.insert(2, sp1[0])
                        op1.insert(3, sp1[1])
                        op1.insert(8, op1[1])
                        op1[4] = "Genome  (-)"
                        new_raw_list1.append(op1)

        elif check_len == 2:
            loca21 = [[500, 750], [820, 1120]]
            if loca21 == []:
                new_raw_list1 = []

            else:
                if list3 != []:
                    for mp, np in zip(list3, loca21):
                        mp.insert(2, np[0])
                        mp.insert(3, np[1])
                        mp.insert(8, mp[1])
                        mp[4] = "Genome  (+)"

                        new_raw_list1.append(mp)
                if list4 != []:
                    for mp1, np1 in zip(list4, loca21):
                        mp1.insert(2, np1[0])
                        mp1.insert(3, np1[1])
                        mp1.insert(8, mp1[1])
                        mp1[4] = "Genome  (-)"
                        new_raw_list1.append(mp1)
        else:
            new_raw_list0 = []
            new_raw_list1 = []
            new_raw_list2 = []

    if new_raw_list0 == [] and new_raw_list1 == [] and new_raw_list2 == []:
        error4 = "<h4>Error</h4><p>The rps12 gene found in your GenBank file does match the two-exon model or three-exon model. As a result, the rps12 gene structure map are not generated.</p><p>This does not necessary mean that the rps12 gene structure is incorrect in your GenBank file</p> "
    else:
        error4 = ""

    record1 = ""
    record2 = ""

    if error1 != "":
        record1 = error1

    elif error1 == "" and error2 != "":
        record1 = error2
    else:
        record1 = ""

    if error1 == "" and error3 != "":
        record2 = error3
    elif error1 == "" and error3 == "" and error4 != "":
        record2 = error4
    else:
        record2 = ""

    if record1 != "" or record2 != "":
        with open("error1.txt", "w") as f:
            f.write(record1)
            f.write(record2)
            print("The error has been output to a file！")

    # The seventh part: save data to CSV file
    folder1 = os.path.exists(file1)
    folder2 = os.path.exists(file2)

    if not folder1 and not folder2:

        os.makedirs(file1)
        print("---  New folder ", file1, " has been generated  ---")

        os.makedirs(file2)
        print("---  New folder ", file2, " has been generated  ---")

        ab_path1 = str(get_path('./' + file1))
        ab_path2 = str(get_path('./' + file2))
        if new_list11 != [] and new_list23 != []:
            name1 = ab_path1 + "/cis_splicing_gene.csv"
            name2 = ab_path1 + "/cis_splicing_subgene.csv"

            csvfile1 = open(name1, 'w', newline='', encoding="utf-8")
            writer1 = csv.writer(csvfile1)
            writer1.writerow(["Gene", "start", "end", "strand", "direction", "Gene_label"])
            writer1.writerows(new_list11)
            csvfile1.close()

            csvfile2 = open(name2, 'w', newline='', encoding="utf-8")
            writer2 = csv.writer(csvfile2)
            writer2.writerow(
                ["Gene", "start", "end", "strand", "direction", "Subgene", "from", "to", "label_from", "label_to",
                 "Gene_label"])
            writer2.writerows(new_list23)
            csvfile2.close()

        if new_raw_list0 != [] and new_raw_list1 != [] and new_raw_list2 != []:
            name3 = ab_path2 + "/trans_splicing_gene.csv"
            name4 = ab_path2 + "/trans_splicing_subgene.csv"

            csvfile3 = open(name3, 'w', newline='', encoding="utf-8")
            writer3 = csv.writer(csvfile3)
            writer3.writerow(
                ["Gene", "Exon", "start", "end", "strand", "direction", "label_from", "label_to", "label_label"])
            writer3.writerows(new_raw_list0)
            writer3.writerows(new_raw_list1)
            writer3.writerows(new_raw_list2)
            csvfile3.close()

            csvfile4 = open(name4, 'w', newline='', encoding="utf-8")
            writer4 = csv.writer(csvfile4)
            writer4.writerow(
                ["Gene", "Exon", "start", "end", "strand", "direction", "label_from", "label_to", "label_label"])
            writer4.writerows(new_raw_list0)
            writer4.writerows(new_raw_list2)
            csvfile4.close()

        elif new_raw_list0 == [] and new_raw_list1 != [] and new_raw_list2 == []:
            namex = ab_path2 + "/trans_splicing_gene.csv"
            namey = ab_path2 + "/trans_splicing_subgene.csv"
            csvfile7 = open(namex, 'w', newline='', encoding="utf-8")
            writer7 = csv.writer(csvfile7)
            writer7.writerow(
                ["Gene", "Exon", "start", "end", "strand", "direction", "label_from", "label_to", "label_label"])
            writer7.writerows(new_raw_list1)
            csvfile7.close()

            csvfile8 = open(namey, 'w', newline='', encoding="utf-8")
            writer8 = csv.writer(csvfile8)
            writer8.writerow(
                ["Gene", "Exon", "start", "end", "strand", "direction", "label_from", "label_to", "label_label"])
            writer8.writerows(new_raw_list1)
            csvfile8.close()
        else:
            print("The structure of rps12 gene does not fit our system")

        print("")
        print("---  All files have been saved to folders!  ---")

    else:
        print(
            "---  New folder ", file1, " already exists here. Please delete it or them first!  ---")
except:
    error0 = "<h4>Error!</h4><p>The corresponding maps cannot be generated. This may be due to the following reasons.</p> <p>1. The uploaded file is not a valid GenBank file;<br>2. The system cannot process the uploaded GenBank file.</p><p>Please check the format and content of the GenBank file</p> "

if error0 != "":
    with open("error1.txt", "w") as f:
        f.write(error0)
        print("The error has been output to a file！")
