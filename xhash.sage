load("drsolve_sage_interface.sage")


def xhash_elimination(dixon_path="./drsolve", live_output=True):
    set_dixon_path(dixon_path)

    field = 4611686018427388039
    ideal = [
        "x1^3 = -7*x0^3 + 8*x0^2 + 8*x0 + 4",
        "x2^3 = -x0^3 + 3*x0^2 - 3*x0 - 5",
        "x3^3 = 8*x0^3 + 5*x0^2 - 8*x0 - 7",
        "x4 = -3*x1^2*x2 + 4*x1^2*x3 - 4*x1^2 + 5*x1*x2^2 + 5*x1*x2*x3 - 8*x1*x2 + x1*x3^2 + 2*x1*x3 - 5*x1 + 5*x2^3 + 6*x2^2*x3 + 6*x2^2 - 7*x2*x3^2 - 5*x2*x3 + 7*x2 - 7*x3^3 + 5*x3^2 - 7*x3 + 7",
        "x5 = -x1^3 - 3*x1^2*x2 - 2*x1^2*x3 + 7*x1^2 - 4*x1*x2^2 + x1*x2*x3 - 6*x1*x2 - 3*x1*x3^2 + 4*x1*x3 + 8*x1 - 3*x2^3 - x2^2*x3 - 2*x2^2 + 7*x2*x3^2 - x2*x3 + 8*x2 + 7*x3^3 + 2*x3^2 - 2*x3",
        "x6 = -7*x1^3 - 7*x1^2*x2 - 2*x1^2*x3 - 7*x1^2 - 2*x1*x2^2 - 3*x1*x2 - x1*x3^2 - 7*x1*x3 - 4*x1 + 3*x2^3 - 3*x2^2*x3 - 4*x2^2 + x2*x3^2 + 3*x2*x3 - 8*x2 + 5*x3^3 - 5*x3 + 5",
    ]

    f0 = "-7*x0^3 + 8*x0^2 + 8*x0 - x1^3 + 4"
    f1 = "-x0^3 + 3*x0^2 - 3*x0 - x2^3 - 5"
    f2 = "8*x0^3 + 5*x0^2 - 8*x0 - x3^3 - 7"
    f3 = "-3*x1^2*x2 + 4*x1^2*x3 - 4*x1^2 + 5*x1*x2^2 + 5*x1*x2*x3 - 8*x1*x2 + x1*x3^2 + 2*x1*x3 - 5*x1 + 5*x2^3 + 6*x2^2*x3 + 6*x2^2 - 7*x2*x3^2 - 5*x2*x3 + 7*x2 - 7*x3^3 + 5*x3^2 - 7*x3 - x4 + 7"
    f4 = "-x1^3 - 3*x1^2*x2 - 2*x1^2*x3 + 7*x1^2 - 4*x1*x2^2 + x1*x2*x3 - 6*x1*x2 - 3*x1*x3^2 + 4*x1*x3 + 8*x1 - 3*x2^3 - x2^2*x3 - 2*x2^2 + 7*x2*x3^2 - x2*x3 + 8*x2 + 7*x3^3 + 2*x3^2 - 2*x3 - x5"
    f5 = "-7*x1^3 - 7*x1^2*x2 - 2*x1^2*x3 - 7*x1^2 - 2*x1*x2^2 - 3*x1*x2 - x1*x3^2 - 7*x1*x3 - 4*x1 + 3*x2^3 - 3*x2^2*x3 - 4*x2^2 + x2*x3^2 + 3*x2*x3 - 8*x2 + 5*x3^3 - 5*x3 - x6 + 5"
    f6 = "-7*x4^3 - 7*x4^2*x5 + 2*x4^2*x6 - 7*x4^2 + 7*x4*x5^2 + 6*x4*x5*x6 + 2*x4*x5 + 8*x4*x6^2 + 5*x4*x6 - x4 + 8*x5^3 - 7*x5^2*x6 - 8*x5^2 - 8*x5*x6^2 - 5*x5*x6 - x6^3 + x6^2 - 7*x6 - 6"

    common = dict(
        field_size=field,
        live_output=live_output,
        verbosity=1,
    )

    print("=== xhash elimination using Sage interface ===")
    r1 = DixonIdeal([f6, f5], ideal, ["x6"], **common)
    #print("r1 = %s" % r1)
    r2 = DixonIdeal([r1, f4], ideal, ["x5"], **common)
    #print("r2 = %s" % r2)
    r3 = DixonIdeal([r2, f3], ideal, ["x4"], **common)
    #print("r3 = %s" % r3)
    r4 = DixonIdeal([r3, f2], ideal, ["x3"], **common)
    #print("r4 = %s" % r4)
    r5 = DixonIdeal([r4, f1], ideal, ["x2"], **common)
    #print("r5 = %s" % r5)
    r6 = DixonIdeal([r5, f0], ideal, ["x1"], **common)
    #print("r6 = %s" % r6)
    return r6


if __name__ == "__main__":
    xhash_elimination()
