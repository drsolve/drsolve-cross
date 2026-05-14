load("drsolve_sage_interface.sage")


def vision_2r_attack(dixon_path="./drsolve", live_output=True):
    set_dixon_path(dixon_path)

    field = "2^8: z8^8 + z8^4 + z8^3 + z8^2 + 1"

    p0 = "(z8^7 + z8^6 + z8^5 + 1)*x1^4*x3^4 + (z8^4 + z8^3 + z8^2 + z8 + 1)*x2^4*x3^4 + (z8^6 + 1)*x1^4*x3^2 + (z8^7 + z8^6 + z8^5 + z8^3 + z8^2 + z8 + 1)*x2^4*x3^2 + (z8^6 + 1)*x1^2*x3^4 + (z8^7 + z8^6 + z8^5 + z8^3 + z8^2 + z8 + 1)*x2^2*x3^4 + (z8^7 + z8^4 + z8^3 + z8^2 + 1)*x1^4*x3 + (z8^6 + z8^4 + z8^3 + z8^2 + 1)*x2^4*x3 + (z8^7 + z8^4 + z8^3 + z8^2 + 1)*x1*x3^4 + (z8^6 + z8^4 + z8^3 + z8^2 + 1)*x2*x3^4 + (z8^7 + z8^4 + z8^2)*x1^2*x3^2 + (z8^7 + z8^6 + z8^4 + z8^3 + z8^2 + z8)*x2^2*x3^2 + (z8^4 + z8^3 + 1)*x3^4 + (z8^7 + z8^6 + z8^4 + z8^3 + z8)*x1^2*x3 + (z8^7 + z8^5 + z8^4 + z8^2 + z8 + 1)*x2^2*x3 + (z8^7 + z8^6 + z8^4 + z8^3 + z8)*x1*x3^2 + (z8^7 + z8^5 + z8^4 + z8^2 + z8 + 1)*x2*x3^2 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^2)*x1*x3 + (z8^7 + z8^3 + z8^2 + z8)*x2*x3 + (z8^6 + z8^5 + z8^3 + z8)*x3^2 + (z8^6 + z8^3 + z8^2 + z8 + 1)*x3 + 1"
    p1 = "(z8^5 + z8^4 + z8^3 + z8^2 + z8)*x1^4*x4^4 + (z8^7 + z8^6)*x2^4*x4^4 + (z8^7 + z8^6 + z8 + 1)*x1^4*x4^2 + (z8^6 + z8^5 + z8^3 + z8^2 + 1)*x2^4*x4^2 + (z8^7 + z8^6 + z8 + 1)*x1^2*x4^4 + (z8^6 + z8^5 + z8^3 + z8^2 + 1)*x2^2*x4^4 + (z8^7 + z8^5 + z8^4 + z8^3 + z8)*x1^4*x4 + (z8^6 + z8^5 + z8^4 + z8^3 + z8)*x2^4*x4 + (z8^7 + z8^5 + z8^4 + z8^3 + z8)*x1*x4^4 + (z8^6 + z8^5 + z8^4 + z8^3 + z8)*x2*x4^4 + (z8^7 + z8^5 + 1)*x1^2*x4^2 + (z8^7 + z8^6 + z8^5 + z8^3 + z8 + 1)*x2^2*x4^2 + (z8^6 + z8^5 + z8 + 1)*x4^4 + (z8^6 + z8^5 + z8^4 + z8 + 1)*x1^2*x4 + (z8^4 + z8^3 + z8^2 + z8)*x2^2*x4 + (z8^6 + z8^5 + z8^4 + z8 + 1)*x1*x4^2 + (z8^4 + z8^3 + z8^2 + z8)*x2*x4^2 + x1*x4 + (z8^6 + z8^5 + z8^4 + z8^3 + z8 + 1)*x2*x4 + (z8^6 + z8^5 + z8^4 + z8^2)*x4^2 + (z8^5 + z8^4 + z8^2)*x4 + 1"
    p2 = "z8*x0*x1 + z8*x1 + 1"
    p3 = "(z8^2 + z8)*x0*x2 + (z8^7 + z8^5 + z8^3 + z8 + 1)*x2 + 1"
    p4 = "z8*x3*x5 + (z8 + 1)*x4*x5 + (z8^5 + z8^4 + z8^2)*x5 + 1"
    p5 = "(z8^2 + z8)*x3*x6 + (z8^2 + z8 + 1)*x4*x6 + (z8^6 + z8^3 + z8^2 + z8)*x6 + 1"
    p6 = "(z8^7 + z8^6 + 1)*x5^4 + (z8^5 + z8^3 + z8^2 + z8 + 1)*x6^4 + (z8^5 + z8^3 + 1)*x5^2 + (z8^7 + z8^5 + z8^4 + z8 + 1)*x6^2 + (z8^6 + z8^5 + z8^4 + z8^3 + z8^2 + 1)*x5 + (z8^7 + z8^6 + z8^3 + z8^2 + 1)*x6 + (z8^6 + z8^4 + z8)"

    common = dict(
        field_size=field,
        live_output=live_output,
        verbosity=1,
        time=True,
    )
    subres = dict(common, resultant_method="subres")
    dixon = dict(common, resultant_method="dixon")

    print("=== Vision Attack using Sage interface ===")
    r1 = DixonRes([p2, p3], ["x0"], **subres)
    #print("r1 = %s" % r1)
    if r1 is None:
        raise RuntimeError("vision r1 failed")

    r2 = DixonRes([r1, p0, p1], ["x1", "x2"], **dixon)
    #print("r2 = %s" % r2)
    if r2 is None:
        raise RuntimeError("vision r2 failed")

    s1 = DixonRes([p6, p5], ["x6"], **subres)
    #print("s1 = %s" % s1)
    if s1 is None:
        raise RuntimeError("vision s1 failed")

    s2 = DixonRes([s1, p4], ["x5"], **subres)
    #print("s2 = %s" % s2)
    if s2 is None:
        raise RuntimeError("vision s2 failed")

    d = DixonRes([r2, s2], ["x3"], **subres)
    #print("d = %s" % d)
    return d


if __name__ == "__main__":
    vision_2r_attack()
