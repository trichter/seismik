# (C) 2018-2020, Tom Eulenfeld, MIT license

import tempfile

def parse(bla):
    if ',' in bla:
        return sum((parse(b) for b in bla.split(',')), tuple())
    if '-' in bla:
        a,b = bla.split('-')
        return tuple(range(int(a), int(b) + 1))
    return (int(bla),)


def read_shotpoints(fname, verbose=False):
    """Load mapping shotpoint: receiver"""
    with open(fname) as f:
        content = f.read()
    sps = {}
    sp = 0
    spacing=None
    for line in content.splitlines():
        if line.startswith('#'):
            continue
        elif '#' in line:
            line = line.split('#')[0].strip()
        spn, term = line.split()
        spn = int(spn)
        while spn > sp + 1:
            sp += 1
            sps[sp] = sps[sp-1] + spacing
        if term[0] == '!':
            spacing = int(term[1:])
        elif term[0] in '+-':
            sps[spn] = sps[sp] + spacing + int(term)
            sp = spn
        elif term == 'X':
            sps[spn] = sps[sp] + spacing
        else:
            # term[0] -> term
            sps[spn] = int(term)
            sp = spn
    if verbose:
        print('\nread shot points - at geofons:')
        print(sps)
    return sps


def read_filenumbers(fname, verbose=False):
    """Load mapping shotpoint: file numbers"""
    with open(fname) as f:
        content = f.read()
    stack = {}
    sp = 0
    default = None
    last_file = 0
    for line in content.splitlines():
        if line.startswith('#'):
            continue
        elif '#' in line:
            line = line.split('#')[0]
        spn, term = line.split()
        spn = int(spn)
        while spn > sp + 1:
            sp += 1
            f1 = last_file + 1
            last_file += default
            stack[sp] = tuple(range(f1, last_file + 1))
        if term[0] == '!':
            default = int(term[1:])
        elif term == 'X':
                f1 = last_file + 1
                f2 = last_file + 1 + default
                stack[spn] = tuple(range(f1, f2))
        else:
            sp += 1
            stack[spn] = parse(term)
            last_file = stack[spn][-1]
    if verbose:
        print('\nread filenumbers for stack')
        print(stack)
    return stack


def read_spreads(fname, verbose=False):
    with open(fname) as f:
        content = f.read()
    spreads = []
    for line in content.splitlines():
        if line.startswith('#'):
            continue
        elif '#' in line:
            line = line.split('#')[0]
        spread, recs, srcs, *special = line.split()
        srcs = parse(srcs)
        recs = parse(recs)
        spread = {'shotpoints': srcs, 'channels': recs}
        for s in special:
            command, r = s.split(':')
            r = parse(r)
            assert command in ('order', 'polarity', 'mute')
            if command == 'order':
                r = (r[0], r[-1])
                try:
                    spread['order'].append(r)
                except KeyError:
                    spread['order'] = [r]
            else:
                spread[command] = r
        spreads.append(spread)
        #for src in parse(srcs):
        #    shot_config[src] = recs
    if verbose:
        print('\nread spreads file: shots, receivers')
        print(spreads)
    return spreads


def test_read_shotpoints():
    test_case = """# test
    1 3
    2 !2
    4 -1
    6 +1
    8 !1
    10 20
    12 X"""
    with tempfile.NamedTemporaryFile('w') as tempf:
        tempf.write(test_case)
        tempf.seek(0)
        read_shotpoints(tempf.name, verbose=True)


def test_read_filenumbers():
    test_case = """# test
    1 3
    2 !2
    4 7-8,10
    6 16-20
    8 !3
    10 X"""
    with tempfile.NamedTemporaryFile('w') as tempf:
        tempf.write(test_case)
        tempf.seek(0)
        read_filenumbers(tempf.name, verbose=True)


def test_read_spreads():
    test_case = """# spread receivers shots
    1      1-12,14-24  1-13"""
    with tempfile.NamedTemporaryFile('w') as tempf:
        tempf.write(test_case)
        tempf.seek(0)
        read_spreads(tempf.name, verbose=True)


if __name__ == '__main__':
    test_read_shotpoints()
    test_read_filenumbers()
    test_read_spreads()
