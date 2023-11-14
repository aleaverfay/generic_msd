import blargs
import re


if __name__ == "__main__":
    with blargs.Parser(locals()) as p:
        p.str("logfile").required()
        p.int("output_index").required()
        p.multiword("requested_vals").cast(lambda x: x.split())
        p.multiword("energies_by_match").cast(lambda x: x.split())

    start_re = re.compile("apps.public.design.mpi_msd: \(0\) Top set")
    outputs_re = re.compile("apps.public.design.mpi_msd: \(0\) Writing structure")
    subexp_re = re.compile(
        "protocols.pack_daemon.DynamicAggregateFunction: \(0\) sub expression"
    )

    if energies_by_match is None:
        energies_by_match = []

    lines = open(logfile).readlines()
    outputs = []
    num = -1
    for line in lines:
        if start_re.match(line):
            if num != -1:
                outputs.append(out)
            out = {}
            out["output_tags"] = []
            num = line.split()[4][1:-1]
        elif outputs_re.match(line):
            output_file = line.split()[4]
            out["output_tags"].append(output_file)
            out[output_file] = line.split()[7]
        elif subexp_re.match(line):
            cols = line.split()
            out[cols[4]] = cols[7]
    outputs.append(out)

    out = outputs[output_index - 1]
    for partial_name in energies_by_match:
        for tag in out["output_tags"]:
            if tag.find(partial_name) != -1:
                print(out[tag], end=" ")
    for key in requested_vals:
        if key not in out:
            print("undefined", end=" ")
        else:
            print(out[key], end=" ")
    print()
