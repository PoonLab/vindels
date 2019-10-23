gp120 = open("/home/jpalmer/PycharmProjects/hiv-withinhost/gp120.txt",'rU')
ntref = ""
for line in gp120:
    line = line.strip("\n")
    ntref += line
c2 = ntref[588:885]
v3 = ntref[885:993]
c3 = ntref[993:1152]

print(c2)
print(v3)
print(c3)

segment = ntref[588:1152]

outfile = open("gp120_c2c3.txt", "w+")
outfile.write(">C2C3_588-1152\n" + segment + "\n" + ">C2_588-885\n" + c2 + "\n"  + ">V3_885-993\n" + v3 + "\n"  + ">C3_993-1152\n" + c3 + "\n")

v_regions = [(390, 468), (468, 588), (885, 993), (1152, 1254), (1377, 1410)]
c_regions = [(0,390),  (588,885) , (993, 1152), (1254, 1377), (1410, 1532)]
