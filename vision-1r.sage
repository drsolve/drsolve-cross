load("drsolve_sage_interface.sage")


def vision_1r_attack(dixon_path="./drsolve", live_output=True):
    set_dixon_path(dixon_path)

    field = "2^8: z8^8 + z8^4 + z8^3 + z8^2 + 1"

    p0 = "z8*x0*x1 + (z8^7 + z8^5 + z8^3 + z8)*x1 + 1"
    p1 = "(z8^2 + z8)*x0*x2 + (z8^6 + z8^5 + z8^4 + z8^3 + z8 + 1)*x2 + 1"
    p2 = "(z8^6 + z8^5 + z8^3)*x1^4 + (z8^6 + z8^4 + z8^3 + z8^2)*x2^4 + (z8^6 + z8^4 + z8^3 + z8^2)*x1^2 + (z8^6 + z8^5 + z8^4 + z8)*x2^2 + (z8^7 + z8^6 + z8^4)*x1 + (z8^7 + z8^5 + z8^4 + z8^3)*x2 + (z8^5 + z8^4 + z8^3 + 1)"

    common = dict(
        field_size=field,
        live_output=live_output,
        verbosity=1,
        time=True,
    )
    subres = dict(common, resultant_method="subres")

    print("=== Vision Attack using Sage interface ===")
    r1 = DixonRes([p0, p1], ["x0"], **subres)
    if r1 is None:
        raise RuntimeError("vision r1 failed")

    d = DixonRes([r1, p2], ["x1"], **subres)
    return d


if __name__ == "__main__":
    vision_1r_attack()
