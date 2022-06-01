IN_INIT = "bondInfoIni.csv"
OUT_INIT = "bondInfoIniFixed.csv"
IN_FIN = "bondInfoFin.csv"
OUT_FIN = "bondInfoFinFixed.csv"

with open(IN_INIT, 'r') as infile:
    with open(OUT_INIT, 'w') as outfile:
        for num, line in enumerate(infile):
            if(num%2 == 0): ##0, 2, 4, 6, 8, ...
                splEven = line.split()
                del splEven[0]
                del splEven[0]
                del splEven[0]
            else: #1, 3, 5, 7, ...
                splOdd = line.split(',')
                wr = splEven + splOdd

                wr_ = ""
                for w in wr:
                   wr_ += str(w) + ','

                outfile.write(wr_[:-1])

        outfile.close()
    infile.close()

with open(IN_FIN, 'r') as infile:
    with open(OUT_FIN, 'w') as outfile:
        for num, line in enumerate(infile):

            spl = line.split(',')
            spl0 = spl.pop(0).split()
            del spl0[0]
            del spl0[0]
            del spl0[0]
            wr = spl0 + spl

            wr_ = ""
            for w in wr:
                wr_ += str(w) + ','

            outfile.write(wr_[:-1])

        outfile.close()
    infile.close()
