# Error bit definitions and descriptions
const ERROR_BIT_DEFINITIONS = [
    (0, 1, "reference array pixels"),
    (1, 2, "reference pixels"),
    (2, 4, "bad reference pixels"),
    (3, 8, "pixels not dark corrected"),
    (4, 16, "pixels with negative dark current"),
    (5, 32, "pixels with large dark current"),
    (6, 64, "flat response too low"),
    (7, 128, "reads dropped for CR rejection = 1"),
    (8, 256, "reads dropped for CR rejection > 1"),
    (9, 512, "bad linear SUTR chi2"),
    (10, 1024, "failed 1D extraction"),
    (11, 2048, "no nearby good pixels in 1D extraction"),
    (12, 4096, "neff>10 in 1D extraction"),
    (13, 8192, "pixel partially saturated"),
    (14, 16384, "pixel fully saturated")
]

# bad_dark_pix_bits = 2^2 + 2^4 #+ 2^5; temporarily remove 2^5 from badlist for now
bad_dark_pix_bits = 2^1 + 2^2 + 2^4
bad_flat_pix_bits = 2^6;
# most multiread CR detections are bad for other reasons
bad_cr_pix_bits = 2^7 + 2^8; # could probably drop 2^7 at least in the future (happily correct 1 read CRs)
bad_chi2_pix_bits = 2^9;

# flags for 1d flux extraction
bad_1d_failed_extract = 2^10;
bad_1d_no_good_pix = 2^11;
bad_1d_neff = 2^12;

bad_pix_bits = bad_dark_pix_bits + bad_flat_pix_bits + bad_cr_pix_bits + bad_chi2_pix_bits +
               bad_1d_failed_extract + bad_1d_no_good_pix + bad_1d_neff;

# Pretty-print error bits
function print_error_bits(error_value::Integer)
    if error_value == 0
        println("No error bits set (value: 0)")
        return
    end

    println("Error bits set for value $error_value:")
    for (bit_num, bit_value, description) in ERROR_BIT_DEFINITIONS
        if (error_value & bit_value) != 0
            println("  Bit $bit_num ($bit_value): $description")
        end
    end
end

# Generate markdown table for error bit definitions
function error_bits_markdown_table()
    io = IOBuffer()

    # Header with alignment specification
    println(io, "| Bit Number | Bit Value | Description |")
    println(io, "|------------|-----------|:------------|")

    # Rows
    for (bit_num, bit_value, description) in ERROR_BIT_DEFINITIONS
        println(io, "| $bit_num | $bit_value | $description |")
    end

    return String(take!(io))
end
