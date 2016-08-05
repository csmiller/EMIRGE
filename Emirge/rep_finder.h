#include <string.h>

inline int ascii2bits(char c) {
    switch(c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return 4;
    }
}

inline int get_bit(unsigned int* data, unsigned int kmer) {
    unsigned int index = kmer / (sizeof(unsigned int) * 8);
    unsigned int shift = kmer % (sizeof(unsigned int) * 8);
    return data[index] & 1 << shift;
}
inline void set_bit(unsigned int* data, unsigned int kmer) {
    unsigned int index = kmer / (sizeof(unsigned int) * 8);
    unsigned int shift = kmer % (sizeof(unsigned int) * 8);
    data[index] |= 1 << shift;
}
inline void unset_bit(unsigned int* data, unsigned int kmer) {
    unsigned int index = kmer / (sizeof(unsigned int) * 8);
    unsigned int shift = kmer % (sizeof(unsigned int) * 8);
    data[index] &= ~(1 << shift);
}

void find_repeats(unsigned int* kmer_map, int k, int min,
                  const char* sequence, int sequence_len,
                  int* results, int maxresults) {
    unsigned int mask = (1L << (2*k)) -1;
    unsigned int kmer = 0;
    int doubles = 0;
    int results_n = 0;
    int pos;

    for (pos = 0; pos < sequence_len; pos++) {
        int base = ascii2bits(sequence[pos]);
        //if (base > 3) continue;
        kmer = ((kmer << 2) | base) & mask;

        if (get_bit(kmer_map, kmer)) {
            doubles ++;
        } else {
	        int match_length = doubles + k - 1;
            if (match_length >= min) {
	            size_t match_start = pos - match_length;

	            char* first = (char*) memmem(
	                sequence, match_start + match_length - 1,
                    sequence + match_start, match_length);

	            if (first) {
	                results[results_n++] = match_length;
	                results[results_n++] = first - sequence;
	                results[results_n++] = match_start;
	                if (results_n + 4 >= maxresults)
	                    break;
	            }
            }
            doubles = 0;
        }
        set_bit(kmer_map, kmer);
    }
    results[results_n++] = 0;

    kmer = 0;
    for (pos = 0; pos < sequence_len; pos++) {
      int base = ascii2bits(sequence[pos]);
      if (base > 3) continue;
      kmer = ((kmer << 2) | base) & mask;
      unset_bit(kmer_map, kmer);
    }
}
