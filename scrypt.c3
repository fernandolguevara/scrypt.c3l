module scrypt;

import std::hash::sha256;

const MAX_BUFFER_LENGTH = (((ulong)1 << 32) - 1) * 32;
const MAX_BLOCKSIZE_PARALLAL_PRODUCT = (ulong)(1 << 30);

macro uint @rotl(uint x, uint n) @local => (((x) << (n)) | ((x) >> (32 - (n))));

macro uint little_endian_u32_at(char[] b, int o) {
    b[o]; // bounds check
    b[o + 3]; // bounds check
    return ((uint)b[o]) | ((uint)(b[o + 1]) << ((uint)8)) | ((uint)(b[o + 2]) << ((uint)16)) | ((uint)(b[o + 3]) << ((uint)24));
}

macro void little_endian_put_u32_at(char[] b, uint v, int o) {
    b[o]; // bounds check
    b[o + 3]; // bounds check
    b[o] = (char)v;
    b[o + 1] = (char)(v >> (uint)(8));
    b[o + 2] = (char)(v >> (uint)(16));
    b[o + 3] = (char)(v >> (uint)(24));
}

macro ulong little_endian_u64_at(char[] b,int o) {
    b[o]; // bounds check
    b[o + 7]; // bounds check
    return ((ulong)b[o]) | (((ulong)b[o + 1]) << ((ulong)8)) | (((ulong)b[o + 2]) << ((ulong)16)) | (((ulong)b[o + 3]) << ((ulong)24)) | (((ulong)b[o + 4]) << ((ulong)32)) | (((ulong)b[o + 5]) << ((ulong)40)) | (((ulong)b[o + 6]) << ((ulong)48)) | (((ulong)b[o + 7]) << ((ulong)56));
}

// salsa20_8 applies the salsa20/8 core transformation to a block
// of 64 u8 bytes.  The block is modified in place.
fn void salsa20_8(char[] block) {
    uint[16] block_words;
    uint[16] scratch;

    int i;
    for (i= 0; i < 16; ++i) {
        block_words[i] = little_endian_u32_at(block, i * 4);
        scratch[i] = block_words[i];
    }

    for (i = 8;i > 0; i -= 2) {    
        // processing columns
        scratch[4] ^= @rotl(scratch[0] + scratch[12], 7);
        scratch[8] ^= @rotl(scratch[4] + scratch[0], 9);
        scratch[12] ^= @rotl(scratch[8] + scratch[4], 13);
        scratch[0] ^= @rotl(scratch[12] + scratch[8], 18);

        scratch[9] ^= @rotl(scratch[5] + scratch[1], 7);
        scratch[13] ^= @rotl(scratch[9] + scratch[5], 9);
        scratch[1] ^= @rotl(scratch[13] + scratch[9], 13);
        scratch[5] ^= @rotl(scratch[1] + scratch[13], 18);

        scratch[14] ^= @rotl(scratch[10] + scratch[6], 7);
        scratch[2] ^= @rotl(scratch[14] + scratch[10], 9);
        scratch[6] ^= @rotl(scratch[2] + scratch[14], 13);
        scratch[10] ^= @rotl(scratch[6] + scratch[2], 18);

        scratch[3] ^= @rotl(scratch[15] + scratch[11], 7);
        scratch[7] ^= @rotl(scratch[3] + scratch[15], 9);
        scratch[11] ^= @rotl(scratch[7] + scratch[3], 13);
        scratch[15] ^= @rotl(scratch[11] + scratch[7], 18);

        // processing rows
        scratch[1] ^= @rotl(scratch[0] + scratch[3], 7);
        scratch[2] ^= @rotl(scratch[1] + scratch[0], 9);
        scratch[3] ^= @rotl(scratch[2] + scratch[1], 13);
        scratch[0] ^= @rotl(scratch[3] + scratch[2], 18);

        scratch[6] ^= @rotl(scratch[5] + scratch[4], 7);
        scratch[7] ^= @rotl(scratch[6] + scratch[5], 9);
        scratch[4] ^= @rotl(scratch[7] + scratch[6], 13);
        scratch[5] ^= @rotl(scratch[4] + scratch[7], 18);

        scratch[11] ^= @rotl(scratch[10] + scratch[9], 7);
        scratch[8] ^= @rotl(scratch[11] + scratch[10], 9);
        scratch[9] ^= @rotl(scratch[8] + scratch[11], 13);
        scratch[10] ^= @rotl(scratch[9] + scratch[8], 18);

        scratch[12] ^= @rotl(scratch[15] + scratch[14], 7);
        scratch[13] ^= @rotl(scratch[12] + scratch[15], 9);
        scratch[14] ^= @rotl(scratch[13] + scratch[12], 13);
        scratch[15] ^= @rotl(scratch[14] + scratch[13], 18);
    }

    for (i= 0; i < 16; ++i) {
        scratch[i] += block_words[i];
        little_endian_put_u32_at(block, scratch[i], i * 4);
    }
}

fn void blkcpy(char[] dest, char[] src, uint len) {
    for (int i = 0; i < len; ++i) {
        dest[i] = src[i];
    }
}


fn void blkxor(char[] dest, char[] src, uint len) {
    for (int i = 0; i < len; ++i) {
        dest[i] ^= src[i];
    }
}

// block_mix performs the block_mix operation using salsa20_8
//
// The block input must be 128 * r in length.  The temp array
// has to be the same size, 128 * r.  r is a positive integer
// value > 0.  The block is modified in place.
fn void block_mix(char[] block, char[] temp, uint r) {
    char[64] scratch;

    blkcpy(&scratch, block[(((2 * r) - 1) * 64)..], 64);

    uint start, stop, temp_start, temp_stop;
    int i;
    for (i = 0; i < 2 * r; ++i) {
        start = i * 64;
        stop = start + 64;

        blkxor(&scratch, block[start..stop - 1], 64);
        salsa20_8(&scratch);

        blkcpy(temp[start..stop - 1], &scratch, 64);
    }

    for (i = 0; i < r; ++i) {
        start = i * 64;
        stop = start + 64;

        temp_start = (i * 2) * 64;
        temp_stop = temp_start + 64;

        blkcpy(block[start..stop - 1], temp[temp_start..temp_stop - 1], 64);
    }

    for (i = 0; i < r; ++i) {
        start = (i + r) * 64;
        stop = start + 64;

        temp_start = ((i * 2) + 1) * 64;
        temp_stop = temp_start + 64;

        blkcpy(block[start..stop - 1], temp[temp_start..temp_stop - 1], 64);
    }
}

fn void smix(char[] block, uint r, ulong n, char[] v_block, char[] temp_block) {
    blkcpy(temp_block, block, 128 * r);
    uint y_start = 128 * r;

    int i;
    ulong v_start, v_stop;
    for (i = 0; i < n; ++i) {
        v_start = (i * (128u64 * r));
        v_stop = (v_start + (128u64 * r));

        blkcpy(v_block[v_start..v_stop - 1], temp_block, 128 * r);
        block_mix(temp_block,temp_block[y_start..], r);
    }

    for (i = 0; i < n; ++i) {
        ulong j = little_endian_u64_at(temp_block, ((2 * r) - 1) * 64) & (n - 1);

        v_start = j * (128u64 * r);
        v_stop = v_start + (128u64 * r);

        blkxor(temp_block, v_block[v_start..v_stop - 1], 128 * r);
        block_mix(temp_block, temp_block[y_start..], r);
    }

    blkcpy(block, temp_block, 128 * r);
}

fault ScryptError {
    OUTPUT_BUFFER_LENGTH_ERROR,
    BLOCKSIZE_PARALLEL_PRODUCT_ERROR,
    CPU_MEMORY_COST_ERROR
}

// scrypt performs password based key derivation using the scrypt algorithm.
//
// The input parameters are:
//
//     password - a slice of bytes which is the password being used to
//         derive the key.  Don't leak this value to anybody.
//     salt - a slice of bytes used to make it harder to crack the key.
//     n - CPU/Memory cost parameter, must be larger than 0, a power of 2,
//         and less than 2^(128 * r / 8).
//     r - block size parameter.
//     p - parallelization parameter, a positive integer less than or
//         equal to ((2^32-1) * hLen) / MFLen where hLen is 32 and
//         MFlen is 128 * r.
//     dk_len - intended output length in octets of the derived key;
//         a positive integer less than or equal to (2^32 - 1) * hLen
//         where hLen is 32.
//
// Reasonable values for n, r, and p are n = 1024, r = 8, p = 16.
fn char[]! scrypt(char[] password, char[] salt, ulong n, uint r, uint p, ulong dk_len) {
    if (dk_len > MAX_BUFFER_LENGTH) {
        return ScryptError.OUTPUT_BUFFER_LENGTH_ERROR?;
    }

    if ((ulong)(r * p) >= MAX_BLOCKSIZE_PARALLAL_PRODUCT) {
        return ScryptError.BLOCKSIZE_PARALLEL_PRODUCT_ERROR?;
    }

    // the following is a sneaky way to determine if a number is a
    // power of 2.  Also, a value of 0 is not allowed.
    if ((n & (n - 1)) != 0 || n == 0) {
        return ScryptError.CPU_MEMORY_COST_ERROR?;
    }

    char[] b = mem::temp_new_array(char, 128u64 * r * p);
    sha256::pbkdf2(password, salt, 1, b);

    char[] xy = mem::temp_new_array(char, 256u64 * r);
    char[] v = mem::temp_new_array(char, 128u64 * r * n);

    int i;
    for (i = 0; i < p; ++i) {
        smix(b[i * 128 * r..], r, n, v, xy);
    }

    char[] result = mem::temp_new_array(char, 128u64 * r * p);
    sha256::pbkdf2(password, b, 1, result);

    return result[:dk_len];
}
