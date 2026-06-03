"""Verifier for TSPTW output validity."""
import sys, math

def verify(input_file, output_file):
    # Read cities
    cities = {}
    with open(input_file) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 5:
                cid, x, y, o, c = int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4])
                cities[cid] = (x, y, o, c)

    # Read solution
    with open(output_file) as f:
        lines = [l.strip() for l in f if l.strip()]
    
    header = lines[0].split()
    k_reported = int(header[0])
    len_reported = int(header[1])
    comp_reported = int(header[2])
    
    tour_ids = [int(lines[i+1]) for i in range(k_reported)]

    # Check unique valid IDs
    seen = set()
    for cid in tour_ids:
        assert cid in cities, f"Invalid city ID: {cid}"
        assert cid not in seen, f"Duplicate city ID: {cid}"
        seen.add(cid)
    assert len(tour_ids) == k_reported, f"Tour count mismatch: {len(tour_ids)} vs {k_reported}"

    # Distance function
    def dist(a, b):
        ax, ay = cities[a][0], cities[a][1]
        bx, by = cities[b][0], cities[b][1]
        return int(math.sqrt((ax-bx)**2 + (ay-by)**2) + 0.5)

    # Simulate tour
    start = tour_ids[0]
    ct = max(0, cities[start][2])  # effective arrival at start
    assert cities[start][2] <= ct <= cities[start][3], f"Start city {start} TW violation: t={ct}, window=[{cities[start][2]},{cities[start][3]}]"
    
    total_dist = 0
    for i in range(1, k_reported):
        prev = tour_ids[i-1]
        cur = tour_ids[i]
        d = dist(prev, cur)
        total_dist += d
        ct += d
        eff = max(ct, cities[cur][2])
        assert eff <= cities[cur][3], f"City {cur} (pos {i}) TW violation: effective={eff}, close={cities[cur][3]}"
        ct = eff

    # Return to start
    d = dist(tour_ids[-1], start)
    total_dist += d
    ct += d

    # Verify reported values
    assert total_dist == len_reported, f"Length mismatch: computed={total_dist}, reported={len_reported}"
    assert ct == comp_reported, f"Completion mismatch: computed={ct}, reported={comp_reported}"

    print(f"VALID: {k_reported} cities, length={len_reported}, completion={comp_reported}")
    return True

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python verify.py input.txt output.txt")
        sys.exit(1)
    try:
        verify(sys.argv[1], sys.argv[2])
    except AssertionError as e:
        print(f"INVALID: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: {e}")
        sys.exit(1)
