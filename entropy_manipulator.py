'''
Entropy Manipulator: given a series of (purely random) binary files, apply bit-level manipulations
to decrease the entropy by a given amount.

The manipulations are:
- for clear: clear one bit in each symbol
- for set: set one bit in each symbol
- for flip: flip one bit in each symbol (xor operation)
- for alternate: alternately clear and set one bit in each symbol
- for correlation: set one bit in each symbol identical to the bit of the previous symbol at an index chosen by the previous symbol

The bit index rule:
- correlation: index = (previous_symbol % symbol_bits)
- others:     index = (current_symbol  % symbol_bits)

The entropy levels supported are:
- M64: 64 bits of entropy per symbol (8 bytes)
- M32: 32 bits of entropy per symbol (4 bytes)
- M16: 16 bits of entropy per symbol (2 bytes)
- M8:  8 bits of entropy per symbol (1 byte)
- M4:  4 bits of entropy per symbol (nibble)
'''

import os
from os.path import isfile, join, basename

bins_path = "/data/ccesare/analysis_dp/bins/mgb_bad/0H/bins"
bits = ["clear", "set", "flip", "alternate", "correlation"]
entropies = ["M4", "M8", "M16", "M32", "M64"]

run_parallel = True

# Set to true to run test mode to check if the operations work correctly
# run for example with: python script.py > output.log to have it written to a file
run_test = False
test_data = bytearray([0b00000001, 0b00100011, 0b01000101, 0b01100111, 0b10001001, 0b10101011, 0b11001101, 0b11101111])


def apply_entropy_manipulation(data: bytearray, bit: str = "clear", entropy: str = "M8"):
    clear_toggle = True  # for alternate

    def modify(value, mask, prev_value=0):
        """
        Apply the manipulation at the bit indicated by 'mask'.
        For 'correlation', copy the bit from prev_value at that position.
        For alternate, toggle between clear and set each call.
        """
        nonlocal clear_toggle
        pos = mask.bit_length() - 1  # bit index from mask (avoid free-variable capture)

        if bit == "clear" or (bit == "alternate" and clear_toggle):
            if bit == "alternate":
                clear_toggle = not clear_toggle
            return value & ~mask

        elif bit == "set" or (bit == "alternate" and not clear_toggle):
            if bit == "alternate":
                clear_toggle = not clear_toggle
            return value | mask

        elif bit == "flip":
            return value ^ mask

        elif bit == "correlation":
            prev_bit = (prev_value >> pos) & 1
            return (value & ~mask) | (prev_bit << pos)

        else:
            raise ValueError(f"Unsupported bit operation: {bit}")

    # ---------- M64: 8 bytes (up to tail shorter than 8 kept consistent) ----------
    if entropy == "M64":
        step = 8
        prev_combined = 0
        for i in range(0, len(data), step):
            remaining = min(step, len(data) - i)  # for tail
            width_bits = 8 * remaining

            combined = 0
            for j in range(remaining):
                combined |= data[i + j] << (8 * (remaining - j - 1))

            # choose index source
            if bit == "correlation":
                pos = prev_combined % width_bits
            else:
                pos = combined % width_bits

            mask = 1 << pos
            result = modify(combined, mask, prev_combined) & ((1 << width_bits) - 1)

            prev_combined = combined
            for j in range(remaining):
                data[i + j] = (result >> (8 * (remaining - j - 1))) & 0xFF

            if run_test:
                print(f"Index {i}: combined={combined:b}, mask_pos={pos}, mask={mask:b}, result={result:b}")

    # ---------- M32: 4 bytes ----------
    elif entropy == "M32":
        step = 4
        prev_combined = 0
        for i in range(0, len(data), step):
            remaining = min(step, len(data) - i)
            width_bits = 8 * remaining

            combined = 0
            for j in range(remaining):
                combined |= data[i + j] << (8 * (remaining - j - 1))

            if bit == "correlation":
                pos = prev_combined % width_bits
            else:
                pos = combined % width_bits

            mask = 1 << pos
            result = modify(combined, mask, prev_combined) & ((1 << width_bits) - 1)

            prev_combined = combined
            for j in range(remaining):
                data[i + j] = (result >> (8 * (remaining - j - 1))) & 0xFF

            if run_test:
                print(f"Index {i}: combined={combined:b}, mask_pos={pos}, mask={mask:b}, result={result:b}")

    # ---------- M16: 2 bytes ----------
    elif entropy == "M16":
        step = 2
        prev_combined = 0
        for i in range(0, len(data), step):
            remaining = min(step, len(data) - i)
            width_bits = 8 * remaining

            combined = 0
            for j in range(remaining):
                combined |= data[i + j] << (8 * (remaining - j - 1))

            if bit == "correlation":
                pos = prev_combined % width_bits
            else:
                pos = combined % width_bits

            mask = 1 << pos
            result = modify(combined, mask, prev_combined) & ((1 << width_bits) - 1)

            prev_combined = combined
            for j in range(remaining):
                data[i + j] = (result >> (8 * (remaining - j - 1))) & 0xFF

            if run_test:
                print(f"Index {i}: combined={combined:b}, mask_pos={pos}, mask={mask:b}, result={result:b}")

    # ---------- M8: 1 byte ----------
    elif entropy == "M8":
        prev_byte = 0
        for i in range(len(data)):
            original = data[i]
            if bit == "correlation":
                pos = prev_byte % 8
            else:
                pos = original % 8

            mask = 1 << pos
            result = modify(original, mask, prev_byte) & 0xFF

            if run_test:
                print(f"Index {i}: data={original:08b}, mask_pos={pos}, mask={mask:08b}, result={result:08b}")

            prev_byte = original
            data[i] = result

    # ---------- M4: nibbles (treat stream as 4-bit symbols in order) ----------
    elif entropy == "M4":
        prev_nib = 0  # previous 4-bit symbol across the nibble stream
        for i in range(len(data)):
            original = data[i]
            hi = (original >> 4) & 0xF
            lo = original & 0xF

            # High nibble processing
            if bit == "correlation":
                pos_prev_hi = prev_nib % 4
            else:
                pos_prev_hi = hi % 4
            mask_hi = (1 << (4 + pos_prev_hi)) & 0xFF  # bit positions 4..7
            tmp = modify(original, mask_hi, prev_nib) & 0xFF

            # Update prev_nib to ORIGINAL high nibble (previous symbol for the low nibble)
            prev_nib = hi

            # Low nibble processing
            if bit == "correlation":
                pos_prev_lo = prev_nib % 4  # prev_nib is original hi at this point
            else:
                pos_prev_lo = lo % 4
            mask_lo = (1 << pos_prev_lo) & 0xFF        # bit positions 0..3
            result = modify(tmp, mask_lo, prev_nib) & 0xFF

            # Update prev_nib to ORIGINAL low nibble (previous for next byte's high nibble)
            prev_nib = lo

            data[i] = result

            if run_test:
                print(
                    f"Index {i}: data={original:08b}, hi={hi:04b}, lo={lo:04b}, "
                    f"pos_prev_hi={pos_prev_hi}, pos_prev_lo={pos_prev_lo}, "
                    f"mask_hi={mask_hi:08b}, mask_lo={mask_lo:08b}, result={result:08b}"
                )

    else:
        raise ValueError(f"Unsupported entropy level: {entropy}")


from concurrent.futures import ProcessPoolExecutor, as_completed

def _process_entropy(bit, entropy, bins, bins_path):
    # (kept for compatibility; unused by main() when run_parallel=True)
    output_dir = join("./output", bit, entropy, "bins")
    os.makedirs(output_dir, exist_ok=True)
    for bin_file in bins:
        print(f"Processing {bin_file} for {bit} with {entropy}", flush=True)

        in_path = join(bins_path, bin_file)
        with open(in_path, "rb") as f:
            data = bytearray(f.read())

        apply_entropy_manipulation(data, bit, entropy)

        with open(join(output_dir, basename(bin_file)), "wb") as out_f:
            out_f.write(data)

    return f"done {bit}/{entropy}"


def _process_file(bit, entropy, bin_file, bins_path, output_root="./output"):
    output_dir = join(output_root, entropy, bit, "bins")
    os.makedirs(output_dir, exist_ok=True)
    print(f"Processing {bin_file} for {bit} with {entropy}", flush=True)

    in_path = join(bins_path, bin_file)
    with open(in_path, "rb") as f:
        data = bytearray(f.read())

    apply_entropy_manipulation(data, bit, entropy)

    out_path = join(output_dir, basename(bin_file))
    with open(out_path, "wb") as out_f:
        out_f.write(data)

    return f"done {bit}/{entropy}/{bin_file}"


def main():
    global bins_path, bits, entropies, run_parallel
    bins = [f for f in os.listdir(bins_path) if isfile(join(bins_path, f))]

    # parallel version:
    if run_parallel:
        tasks = [(bit, entropy, bin_file) for bit in bits for entropy in entropies for bin_file in bins]
        max_workers = min(len(tasks), os.cpu_count() or 1)
        with ProcessPoolExecutor(max_workers=max_workers) as exe:
            futures = {
                exe.submit(_process_file, bit, entropy, bin_file, bins_path, "./output"): (bit, entropy, bin_file)
                for bit, entropy, bin_file in tasks
            }
            for fut in as_completed(futures):
                bit_done, entropy_done, bin_done = futures[fut]
                try:
                    print(fut.result(), flush=True)
                except Exception as e:
                    print(f"Error processing {bit_done}/{entropy_done}/{bin_done}: {e}", flush=True)
    else:
        # non-parallel version
        for bit in bits:
            for entropy in entropies:
                output_dir = join("./output", entropy, bit, "bins")
                os.makedirs(output_dir, exist_ok=True)
                for bin_file in bins:
                    print(f"Processing {bin_file} for {bit} with {entropy}", flush=True)
                    with open(join(bins_path, bin_file), "rb") as file:
                        data = bytearray(file.read())

                    apply_entropy_manipulation(data, bit, entropy)
                    out_path = join(output_dir, basename(bin_file))
                    with open(out_path, "wb") as out_f:
                        out_f.write(data)


def test_main():
    global bits, entropies, run_parallel, test_data
    for bit in bits:
        for entropy in entropies:
            print(f"Testing {bit} with {entropy}")
            data = bytearray(test_data)
            apply_entropy_manipulation(data, bit, entropy)

            print("Original \\ Result:")
            print([f"0b{b:08b}" for b in test_data])
            print([f"0b{b:08b}" for b in data])
            print("")


if __name__ == "__main__":
    if run_test:
        test_main()
    else:
        main()
