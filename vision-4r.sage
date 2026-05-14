load("drsolve_sage_interface.sage")


def vision_4r_attack(dixon_path="./drsolve", live_output=True):
    set_dixon_path(dixon_path)

    field = "2^8: z8^8 + z8^4 + z8^3 + z8^2 + 1"

    p1 = "(z8^4 + z8 + 1)*x1^4*x3^4 + (z8^7 + z8^4 + z8^2)*x2^4*x3^4 + (z8^7 + z8^5 + z8^3 + z8^2 + z8)*x1^4*x3^2 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^3 + 1)*x2^4*x3^2 + (z8^7 + z8^5 + z8^3 + z8^2 + z8)*x1^2*x3^4 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^3 + 1)*x2^2*x3^4 + (z8^6 + z8^5 + z8^3 + z8^2 + z8)*x1^4*x3 + (z8^6 + z8^4 + z8^3 + 1)*x2^4*x3 + (z8^6 + z8^5 + z8^3 + z8^2 + z8)*x1*x3^4 + (z8^6 + z8^4 + z8^3 + 1)*x2*x3^4 + (z8^7 + z8^6 + z8^3 + z8^2 + z8)*x1^2*x3^2 + (z8^7 + z8^5 + z8^3 + 1)*x2^2*x3^2 + (z8^7 + z8^6 + z8^4 + z8^3 + 1)*x3^4 + (z8^6 + z8^5 + z8^3)*x1^2*x3 + (z8^6 + z8^4 + z8^3 + z8^2)*x2^2*x3 + (z8^6 + z8^5 + z8^3)*x1*x3^2 + (z8^6 + z8^4 + z8^3 + z8^2)*x2*x3^2 + (z8^7 + z8^3)*x1*x3 + (z8^7 + z8^6 + z8^3 + z8^2)*x2*x3 + (z8^7 + z8^4 + z8^3)*x3^2 + (z8^7 + z8^2 + z8 + 1)*x3 + 1"
    p2 = "(z8^5 + z8^4 + z8^2 + 1)*x1^4*x4^4 + (z8^7 + z8^5 + z8^4 + z8)*x2^4*x4^4 + (z8^7 + z8^6 + z8^5 + z8^3 + z8^2 + z8 + 1)*x1^4*x4^2 + (z8^7 + z8^5 + z8^4 + z8^3)*x2^4*x4^2 + (z8^7 + z8^6 + z8^5 + z8^3 + z8^2 + z8 + 1)*x1^2*x4^4 + (z8^7 + z8^5 + z8^4 + z8^3)*x2^2*x4^4 + (z8^7 + z8^5 + z8^4 + z8)*x1^4*x4 + (z8^7 + z8^2 + 1)*x2^4*x4 + (z8^7 + z8^5 + z8^4 + z8)*x1*x4^4 + (z8^7 + z8^2 + 1)*x2*x4^4 + (z8^6 + z8^3 + z8^2 + z8 + 1)*x1^2*x4^2 + (z8^5 + z8^3)*x2^2*x4^2 + (z8^7 + z8^6 + z8^5 + z8^3)*x4^4 + (z8^7 + z8^5 + z8^4 + z8^3)*x1^2*x4 + (z8^7 + z8^3 + z8^2)*x2^2*x4 + (z8^7 + z8^5 + z8^4 + z8^3)*x1*x4^2 + (z8^7 + z8^3 + z8^2)*x2*x4^2 + (z8^7 + z8^2 + 1)*x1*x4 + (z8^7 + z8^6 + 1)*x2*x4 + (z8^7 + z8^6 + z8^4 + z8^3 + z8^2)*x4^2 + (z8^7 + z8^6 + z8^2)*x4 + 1"
    p3 = "(z8^4 + z8 + 1)*x5^4*x7^4 + (z8^7 + z8^4 + z8^2)*x6^4*x7^4 + (z8^7 + z8^5 + z8^3 + z8^2 + z8)*x5^4*x7^2 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^3 + 1)*x6^4*x7^2 + (z8^7 + z8^5 + z8^3 + z8^2 + z8)*x5^2*x7^4 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^3 + 1)*x6^2*x7^4 + (z8^6 + z8^5 + z8^3 + z8^2 + z8)*x5^4*x7 + (z8^6 + z8^4 + z8^3 + 1)*x6^4*x7 + (z8^6 + z8^5 + z8^3 + z8^2 + z8)*x5*x7^4 + (z8^6 + z8^4 + z8^3 + 1)*x6*x7^4 + (z8^7 + z8^6 + z8^3 + z8^2 + z8)*x5^2*x7^2 + (z8^7 + z8^5 + z8^3 + 1)*x6^2*x7^2 + (z8^7 + z8^6 + z8^5 + z8^4 + 1)*x7^4 + (z8^6 + z8^5 + z8^3)*x5^2*x7 + (z8^6 + z8^4 + z8^3 + z8^2)*x6^2*x7 + (z8^6 + z8^5 + z8^3)*x5*x7^2 + (z8^6 + z8^4 + z8^3 + z8^2)*x6*x7^2 + (z8^7 + z8^3)*x5*x7 + (z8^7 + z8^6 + z8^3 + z8^2)*x6*x7 + (z8^7 + z8^6 + z8^5 + z8)*x7^2 + (z8^7 + z8^5 + z8 + 1)*x7 + 1"
    p4 = "(z8^5 + z8^4 + z8^2 + 1)*x5^4*x8^4 + (z8^7 + z8^5 + z8^4 + z8)*x6^4*x8^4 + (z8^7 + z8^6 + z8^5 + z8^3 + z8^2 + z8 + 1)*x5^4*x8^2 + (z8^7 + z8^5 + z8^4 + z8^3)*x6^4*x8^2 + (z8^7 + z8^6 + z8^5 + z8^3 + z8^2 + z8 + 1)*x5^2*x8^4 + (z8^7 + z8^5 + z8^4 + z8^3)*x6^2*x8^4 + (z8^7 + z8^5 + z8^4 + z8)*x5^4*x8 + (z8^7 + z8^2 + 1)*x6^4*x8 + (z8^7 + z8^5 + z8^4 + z8)*x5*x8^4 + (z8^7 + z8^2 + 1)*x6*x8^4 + (z8^6 + z8^3 + z8^2 + z8 + 1)*x5^2*x8^2 + (z8^5 + z8^3)*x6^2*x8^2 + (z8^6 + z8^5 + z8^4 + z8^3 + z8^2 + z8)*x8^4 + (z8^7 + z8^5 + z8^4 + z8^3)*x5^2*x8 + (z8^7 + z8^3 + z8^2)*x6^2*x8 + (z8^7 + z8^5 + z8^4 + z8^3)*x5*x8^2 + (z8^7 + z8^3 + z8^2)*x6*x8^2 + (z8^7 + z8^2 + 1)*x5*x8 + (z8^7 + z8^6 + 1)*x6*x8 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^3 + z8^2 + z8)*x8^2 + (z8^6 + z8^5 + z8^3 + z8 + 1)*x8 + 1"
    p5 = "(z8^4 + z8 + 1)*x9^4*x11^4 + (z8^7 + z8^4 + z8^2)*x10^4*x11^4 + (z8^7 + z8^5 + z8^3 + z8^2 + z8)*x9^4*x11^2 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^3 + 1)*x10^4*x11^2 + (z8^7 + z8^5 + z8^3 + z8^2 + z8)*x9^2*x11^4 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^3 + 1)*x10^2*x11^4 + (z8^6 + z8^5 + z8^3 + z8^2 + z8)*x9^4*x11 + (z8^6 + z8^4 + z8^3 + 1)*x10^4*x11 + (z8^6 + z8^5 + z8^3 + z8^2 + z8)*x9*x11^4 + (z8^6 + z8^4 + z8^3 + 1)*x10*x11^4 + (z8^7 + z8^6 + z8^3 + z8^2 + z8)*x9^2*x11^2 + (z8^7 + z8^5 + z8^3 + 1)*x10^2*x11^2 + (z8^7 + z8^5 + z8 + 1)*x11^4 + (z8^6 + z8^5 + z8^3)*x9^2*x11 + (z8^6 + z8^4 + z8^3 + z8^2)*x10^2*x11 + (z8^6 + z8^5 + z8^3)*x9*x11^2 + (z8^6 + z8^4 + z8^3 + z8^2)*x10*x11^2 + (z8^7 + z8^3)*x9*x11 + (z8^7 + z8^6 + z8^3 + z8^2)*x10*x11 + (z8^7 + z8^6 + z8^3 + z8^2 + 1)*x11^2 + (z8^4 + z8^3 + z8^2 + 1)*x11 + 1"
    p6 = "(z8^5 + z8^4 + z8^2 + 1)*x9^4*x12^4 + (z8^7 + z8^5 + z8^4 + z8)*x10^4*x12^4 + (z8^7 + z8^6 + z8^5 + z8^3 + z8^2 + z8 + 1)*x9^4*x12^2 + (z8^7 + z8^5 + z8^4 + z8^3)*x10^4*x12^2 + (z8^7 + z8^6 + z8^5 + z8^3 + z8^2 + z8 + 1)*x9^2*x12^4 + (z8^7 + z8^5 + z8^4 + z8^3)*x10^2*x12^4 + (z8^7 + z8^5 + z8^4 + z8)*x9^4*x12 + (z8^7 + z8^2 + 1)*x10^4*x12 + (z8^7 + z8^5 + z8^4 + z8)*x9*x12^4 + (z8^7 + z8^2 + 1)*x10*x12^4 + (z8^6 + z8^3 + z8^2 + z8 + 1)*x9^2*x12^2 + (z8^5 + z8^3)*x10^2*x12^2 + (z8^7 + z8^6 + z8^4 + z8^2 + 1)*x12^4 + (z8^7 + z8^5 + z8^4 + z8^3)*x9^2*x12 + (z8^7 + z8^3 + z8^2)*x10^2*x12 + (z8^7 + z8^5 + z8^4 + z8^3)*x9*x12^2 + (z8^7 + z8^3 + z8^2)*x10*x12^2 + (z8^7 + z8^2 + 1)*x9*x12 + (z8^7 + z8^6 + 1)*x10*x12 + (z8^6 + z8^5 + z8^4 + z8^3)*x12^2 + (z8^7 + z8^3 + 1)*x12 + 1"
    p7 = "z8*x0*x1 + (z8^7 + z8^6 + z8^4 + z8^3 + z8 + 1)*x1 + 1"
    p8 = "(z8^2 + z8)*x0*x2 + (z8^7 + z8^6 + z8^4 + z8^3 + z8 + 1)*x2 + 1"
    p9 = "z8*x3*x5 + (z8 + 1)*x4*x5 + (z8^7 + z8^6 + z8^4 + z8^2 + 1)*x5 + 1"
    p10 = "(z8^2 + z8)*x3*x6 + (z8^2 + z8 + 1)*x4*x6 + (z8^7 + z8^3 + z8^2 + z8 + 1)*x6 + 1"
    p11 = "z8*x7*x9 + (z8 + 1)*x8*x9 + (z8^7 + z8^2 + 1)*x9 + 1"
    p12 = "(z8^2 + z8)*x7*x10 + (z8^2 + z8 + 1)*x8*x10 + (z8^7 + z8^5 + z8^3 + z8^2)*x10 + 1"
    p13 = "z8*x11*x13 + (z8 + 1)*x12*x13 + (z8^6 + z8^5 + z8^3)*x13 + 1"
    p14 = "(z8^2 + z8)*x11*x14 + (z8^2 + z8 + 1)*x12*x14 + (z8^7 + z8^6 + z8^3)*x14 + 1"
    p15 = "(z8^7 + z8^5 + z8^3 + 1)*x13^4 + (z8^6 + z8^5 + z8^4 + z8 + 1)*x14^4 + (z8^6 + z8^4 + z8^3 + z8^2 + 1)*x13^2 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^3 + z8^2 + 1)*x14^2 + (z8^4 + z8^2)*x13 + (z8^4 + z8^3 + z8^2 + z8)*x14 + z8^4"

    common = dict(
        field_size=field,
        live_output=live_output,
        verbosity=1,
        time=True,
    )
    subres = dict(common, resultant_method="subres")
    dixon = dict(common, resultant_method="dixon")

    print("=== Vision Attack using Sage interface ===")
    r1 = DixonRes([p7, p8], ["x0"], **subres)
    if r1 is None:
        raise RuntimeError("vision r1 failed")

    r2 = DixonRes([r1, p1, p2], ["x1", "x2"], **dixon)
    if r2 is None:
        raise RuntimeError("vision r2 failed")

    r3 = DixonRes([r2, p9, p10], ["x3", "x4"], **dixon)
    if r3 is None:
        raise RuntimeError("vision r3 failed")

    r4 = DixonRes([r3, p3, p4], ["x5", "x6"], **dixon)
    if r4 is None:
        raise RuntimeError("vision r4 failed")

    s1 = DixonRes([p15, p14], ["x14"], **subres)
    if s1 is None:
        raise RuntimeError("vision s1 failed")

    s2 = DixonRes([s1, p13], ["x13"], **subres)
    if s2 is None:
        raise RuntimeError("vision s2 failed")

    s3 = DixonRes([s2, p6], ["x12"], **subres)
    if s3 is None:
        raise RuntimeError("vision s3 failed")

    s4 = DixonRes([s3, p5], ["x11"], **subres)
    if s4 is None:
        raise RuntimeError("vision s4 failed")

    s5 = DixonRes([s4, p12], ["x10"], **subres)
    if s5 is None:
        raise RuntimeError("vision s5 failed")

    s6 = DixonRes([s5, p11], ["x9"], **subres)
    if s6 is None:
        raise RuntimeError("vision s6 failed")

    d = DixonRes([r4, s6], ["x7"], **subres)
    return d


if __name__ == "__main__":
    vision_4r_attack()
