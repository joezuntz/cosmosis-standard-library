import scipy.special

template = """#ifndef BESSELZEROES_H
#define BESSELZEROES_H

const double BesselJ0_Zeros[100000]={{
{j0_zeros}
}};

const double BesselJ4_Zeros[100000]={{
{j4_zeros}
}};

#endif
"""

def main():
    # use scipy.special to calculate the j0 and j4 zeros
    j0_zeros = scipy.special.jn_zeros(0, 100000)
    j4_zeros = scipy.special.jn_zeros(4, 100000)

    # write the template to a file
    with open("besselzero04.h", "w") as f:
        j0_text = ",\n".join(map(str, j0_zeros))
        j4_text = ",\n".join(map(str, j4_zeros))

        f.write(template.format(
            j0_zeros=j0_text,
            j4_zeros=j4_text,
        ))


if __name__ == "__main__":
    main()
