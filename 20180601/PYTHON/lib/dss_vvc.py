import re
from lib.vvc import voltVARControl
def generate_vvc_objects(fn):
    vvarobjects = []
    f = open(fn, "r")

    for l in f:
        if re.findall('(?i)New VVC', l):
            bp = []

            breakpoints = str.split(re.findall('(?i)BP=\s*([^\n\r]*)', l)[0], ',')

            for elem in breakpoints:
                point = re.findall("[0-9.]*",  elem)
                for i in point:
                    if i:
                        bp.append(float(i))
            print(bp)

            phase = int(re.findall('(?i)phase=\s*([0-9]*)', l)[0]) - 1
            print(phase)

            minkvar = float(re.findall('(?i)min_kvar=([-0-9]*)', l)[0])
            print(minkvar)

            maxkvar = float(re.findall('(?i)max_kvar=([-0-9]*)', l)[0])
            print(maxkvar)

            bus = re.findall('(?i)bus=([\w.]*)\s', l)[0]
            bus = re.findall('[\w]*', bus)[0]
            print(bus)
            voltvarobject = voltVARControl(bp, minkvar, maxkvar, bus, phase)
            vvarobjects.append(voltvarobject)
            print("\n --------------")

    for e in vvarobjects:
        print(e)
    return vvarobjects