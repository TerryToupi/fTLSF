#ifndef __FTLSF__
#define __FTLSF__

#define _t(name) ftlsf_##name

typedef unsigned char _t(u8);
typedef signed char _t(i8);
typedef signed char _t(b8);
typedef unsigned short _t(u16);
typedef signed short _t(i16);
typedef signed short _t(b16);
typedef unsigned int _t(u32);
typedef signed int _t(i32);
typedef signed int _t(b32);
typedef unsigned long long _t(u64);
typedef signed long long _t(i64);
typedef signed long long _t(b64);
typedef float _t(f32);
typedef double _t(f64);

#ifdef __cplusplus
#define FTLSF_STATIC_ASSERT(condition, message) static_assert(condition, message)
#else
#define FTLSF_STATIC_ASSERT(condition, message) _Static_assert(condition, message)
#endif

FTLSF_STATIC_ASSERT((_t(u8))-1 == 0xffu, "fTLSF requires 8-bit bytes");
FTLSF_STATIC_ASSERT(sizeof(_t(u16)) == 2 && (_t(u16))-1 == 0xffffu, "fTLSF requires a 16-bit short");
FTLSF_STATIC_ASSERT(sizeof(_t(u32)) == 4 && (_t(u32))-1 == 0xffffffffu, "fTLSF requires a 32-bit int");
FTLSF_STATIC_ASSERT(sizeof(_t(u64)) == 8 && (_t(u64))-1 == 0xffffffffffffffffull, "fTLSF requires a 64-bit long long");

#undef FTLSF_STATIC_ASSERT

enum
{
  FTLSF_NUM_TOP_BINS         = 32ul,
  FTLSF_BINS_PER_LEAF        = 8ul,
  FTLSF_TOP_BINS_INDEX_SHIFT = 3ul,
  FTLSF_LEAF_BINS_INDEX_MASK = 0x7ul,
  FTLSF_NUM_LEAF_BINS        = FTLSF_NUM_TOP_BINS * FTLSF_BINS_PER_LEAF,
};

typedef struct ftlsf_allocation ftlsf_allocation;
struct ftlsf_allocation
{
  _t(u32) offset;
  _t(u32) meta;
};

typedef struct ftlsf_node ftlsf_node;
struct ftlsf_node
{
  _t(u32) data_offset;
  _t(u32) data_size;
  _t(u32) bin_prev;
  _t(u32) bin_next;
  _t(u32) bin_n_prev;
  _t(u32) bin_n_next;
  _t(b32) used;
};

typedef struct ftlsf_allocator ftlsf_allocator;
struct ftlsf_allocator
{
  _t(u32) size;
  _t(u32) max_allocations;
  _t(u32) used_bins_top;
  _t(u32) free_storage;
  _t(u32) free_offset; /* Number of available entries in free_nodes. */

  _t(u8) used_bins[FTLSF_NUM_TOP_BINS];
  _t(u32) bin_indices[FTLSF_NUM_LEAF_BINS];
  _t(u32) *free_nodes;
  ftlsf_node *nodes;
};

ftlsf_allocator *ftlsf_create_allocator(_t(u32) size, _t(u32) max_allocations);
void ftlsf_destroy_allocator(ftlsf_allocator *allocator);
ftlsf_allocation ftlsf_allocate(ftlsf_allocator *allocator, _t(u64) size);
void ftlsf_free(ftlsf_allocator *allocator, ftlsf_allocation allocation);

#undef _t
#endif


/* IMPLEMENTATION */


#ifdef FTLSF_IMPLEMENTATION
#undef FTLSF_IMPLEMENTATION

#define _t(name) ftlsf_##name

#ifdef DEBUG
#include <assert.h>
#define FTLSF_ASSERT(x) assert(x)
#else
#define FTLSF_ASSERT(x)
#endif

#if defined(_MSC_VER) && !defined(FTLSF_NO_INTRINSICS)
#include <intrin.h>
#endif

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static const _t(u32) ftlsf_mantissa_bits  = 3;
static const _t(u32) ftlsf_mantissa_value = 1 << ftlsf_mantissa_bits;
static const _t(u32) ftlsf_mantissa_mask  = ftlsf_mantissa_value - 1;

/* Define FTLSF_NO_INTRINSICS to force the portable bit-scan implementations. */
#if !defined(FTLSF_NO_INTRINSICS)
#if defined(__has_builtin)
#if __has_builtin(__builtin_clz)
#define FTLSF_HAS_BUILTIN_CLZ
#endif
#if __has_builtin(__builtin_ctz)
#define FTLSF_HAS_BUILTIN_CTZ
#endif
#elif defined(__GNUC__)
#define FTLSF_HAS_BUILTIN_CLZ
#define FTLSF_HAS_BUILTIN_CTZ
#endif
#endif

/*
 * Count zero bits before the highest set bit in v, returning 0..31. Requires v != 0.
 * Uses a compiler intrinsic when available; otherwise fills the lower bits and
 * hashes the result into a De Bruijn lookup table. The product wraps to 32 bits.
 */
static inline _t(u32)
ftlsf_lzcnt_nonzero(_t(u32) v)
{
#if defined(_MSC_VER) && !defined(FTLSF_NO_INTRINSICS)
  unsigned long val;
  _BitScanReverse(&val, v);
  return 31 - val;
#elif defined(FTLSF_HAS_BUILTIN_CLZ)
  return __builtin_clz(v);
#else
  static const _t(u8) leading_zeros[32] = {
    31,
    22,
    30,
    21,
    18,
    10,
    29,
    2,
    20,
    17,
    15,
    13,
    9,
    6,
    28,
    1,
    23,
    19,
    11,
    3,
    16,
    14,
    7,
    24,
    12,
    4,
    8,
    25,
    5,
    26,
    27,
    0,
  };

  /* Fill every bit below the highest set bit before hashing. */
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  return leading_zeros[(v * 0x07c4acddu) >> 27];
#endif
}

/*
 * Count zero bits below the lowest set bit in v, returning 0..31. Requires v != 0.
 * Uses a compiler intrinsic when available; otherwise isolates the lowest set bit
 * and hashes it into a De Bruijn lookup table. The product wraps to 32 bits.
 */
static inline _t(u32)
ftlsf_tzcnt_nonzero(_t(u32) v)
{
#if defined(_MSC_VER) && !defined(FTLSF_NO_INTRINSICS)
  unsigned long val;
  _BitScanForward(&val, v);
  return val;
#elif defined(FTLSF_HAS_BUILTIN_CTZ)
  return __builtin_ctz(v);
#else
  static const _t(u8) trailing_zeros[32] = {
    0,
    1,
    28,
    2,
    29,
    14,
    24,
    3,
    30,
    22,
    20,
    15,
    25,
    17,
    4,
    8,
    31,
    27,
    13,
    23,
    21,
    19,
    16,
    7,
    26,
    12,
    18,
    6,
    11,
    5,
    10,
    9,
  };

  /* Isolate the lowest set bit before hashing. */
  _t(u32) lowest_bit = v & (~v + 1u);
  return trailing_zeros[(lowest_bit * 0x077cb531u) >> 27];
#endif
}

#undef FTLSF_HAS_BUILTIN_CLZ
#undef FTLSF_HAS_BUILTIN_CTZ

/*
 * Encode size as an integer exponent/mantissa bin code, rounding up to the next
 * representable size when needed. Sizes 0..7 are exact; larger sizes use three
 * mantissa bits and an implicit leading bit. Used to select allocation bins.
 * At the upper end of the u32 range, the code can represent a size of 2^32.
 */
static inline _t(u32)
ftlsf_u32_to_float_round_up(_t(u32) size)
{
  _t(u32) exponent = 0;
  _t(u32) mantissa = 0;

  if(size < ftlsf_mantissa_value)
  {
    mantissa = size;
  }
  else
  {
    _t(u32) leading_zeros   = ftlsf_lzcnt_nonzero(size);
    _t(u32) highest_set_bit = 31 - leading_zeros;

    _t(u32) mantissa_start_bit = highest_set_bit - ftlsf_mantissa_bits;
    exponent                   = mantissa_start_bit + 1;
    mantissa                   = (size >> mantissa_start_bit) & ftlsf_mantissa_mask;

    _t(u32) low_bits_mask = (1 << mantissa_start_bit) - 1;

    if((size & low_bits_mask) != 0)
      ++mantissa;
  }

  return (exponent << ftlsf_mantissa_bits) + mantissa;
}

/*
 * Encode size as an integer exponent/mantissa bin code, rounding down to the
 * nearest representable size. Sizes 0..7 are exact. Used to classify free blocks
 * so a bin's represented size never exceeds the actual size of its blocks.
 */
static inline _t(u32)
ftlsf_u32_to_float_round_down(_t(u32) size)
{
  _t(u32) exponent = 0;
  _t(u32) mantissa = 0;

  if(size < ftlsf_mantissa_value)
  {
    mantissa = size;
  }
  else
  {
    _t(u32) leading_zeros   = ftlsf_lzcnt_nonzero(size);
    _t(u32) highest_set_bit = 31 - leading_zeros;

    _t(u32) mantissa_start_bit = highest_set_bit - ftlsf_mantissa_bits;
    exponent                   = mantissa_start_bit + 1;
    mantissa                   = (size >> mantissa_start_bit) & ftlsf_mantissa_mask;
  }

  return (exponent << ftlsf_mantissa_bits) + mantissa;
}

/*
 * Decode an integer exponent/mantissa bin code into its represented size.
 * This is the bin's size boundary, which may differ from the original block size.
 * Requires a code whose decoded size fits in u32.
 */
static inline _t(u32)
ftlsf_float_to_u32(_t(u32) float_value)
{
  _t(u32) exponent = float_value >> ftlsf_mantissa_bits;
  _t(u32) mantissa = float_value & ftlsf_mantissa_mask;

  if(exponent == 0)
  {
    return mantissa;
  }
  else
  {
    return (mantissa | ftlsf_mantissa_value) << (exponent - 1);
  }
}

/*
 * Find the lowest set bit at or above start_bit_index, including that position.
 * Returns its index, or UINT32_MAX if no qualifying bit is set.
 * Requires start_bit_index < 32; bit_mask may be zero.
 */
static inline _t(u32)
ftlsf_find_lowest_set_bit_after(_t(u32) bit_mask, _t(u32) start_bit_index)
{
  _t(u32) mask_before_start_index = (UINT32_C(1) << start_bit_index) - 1;
  _t(u32) mask_after_start_index  = ~mask_before_start_index;
  _t(u32) bits_after              = bit_mask & mask_after_start_index;
  if(bits_after == 0)
    return UINT32_MAX;
  return ftlsf_tzcnt_nonzero(bits_after);
}

/*
 * Round n up to the nearest multiple, leaving an already aligned value unchanged.
 * Supports any nonzero multiple, including values that are not powers of two.
 * Requires the rounded result to fit in u64.
 */
static inline _t(u64)
ftlsf_align_up(_t(u64) n, _t(u64) multiple)
{
  return n + (multiple - n % multiple) % multiple;
}

/*
 * Register a free block of size bytes at data_offset using a spare node record.
 * Push it onto the appropriate bin, set occupancy bits, and add to free_storage.
 * Physical-neighbor links start empty and must be connected by the caller.
 * Returns the node index, or UINT32_MAX when no spare node is available.
 */
static inline _t(u32)
ftlsf_insert_node_into_bin(ftlsf_allocator *allocator, _t(u32) size, _t(u32) data_offset)
{
  if(allocator->free_offset == 0)
    return UINT32_MAX;

  _t(u32) bin_index      = ftlsf_u32_to_float_round_down(size);
  _t(u32) top_bin_index  = bin_index >> FTLSF_TOP_BINS_INDEX_SHIFT;
  _t(u32) leaf_bin_index = bin_index & FTLSF_LEAF_BINS_INDEX_MASK;

  if(allocator->bin_indices[bin_index] == UINT32_MAX)
  {
    allocator->used_bins[top_bin_index] |= 1 << leaf_bin_index;
    allocator->used_bins_top |= 1 << top_bin_index;
  }

  _t(u32) top_node_index = allocator->bin_indices[bin_index];
  _t(u32) node_index     = allocator->free_nodes[--allocator->free_offset];

  allocator->nodes[node_index] = (ftlsf_node){
    .data_offset = data_offset,
    .data_size   = size,
    .bin_next    = top_node_index,
    .bin_prev    = UINT32_MAX,
    .bin_n_next  = UINT32_MAX,
    .bin_n_prev  = UINT32_MAX,
    .used        = 0,
  };

  if(top_node_index != UINT32_MAX)
    allocator->nodes[top_node_index].bin_prev = node_index;
  allocator->bin_indices[bin_index] = node_index;

  allocator->free_storage += size;

  return node_index;
}

/*
 * Unlink a free node from its bin, clear occupancy bits for any emptied bins,
 * subtract its size from free_storage, and return its record to the spare stack.
 * Requires node_index to identify a free node currently linked in a bin.
 * The caller is responsible for reconnecting its physical neighbors.
 */
static inline void
ftlsf_remove_node_from_bin(ftlsf_allocator *allocator, _t(u32) node_index)
{
  ftlsf_node *node = &allocator->nodes[node_index];

  if(node->bin_prev != UINT32_MAX)
  {
    allocator->nodes[node->bin_prev].bin_next = node->bin_next;
    if(node->bin_next != UINT32_MAX)
      allocator->nodes[node->bin_next].bin_prev = node->bin_prev;
  }
  else
  {
    _t(u32) bin_index      = ftlsf_u32_to_float_round_down(node->data_size);
    _t(u32) top_bin_index  = bin_index >> FTLSF_TOP_BINS_INDEX_SHIFT;
    _t(u32) leaf_bin_index = bin_index & FTLSF_LEAF_BINS_INDEX_MASK;

    allocator->bin_indices[bin_index] = node->bin_next;
    if(node->bin_next != UINT32_MAX)
      allocator->nodes[node->bin_next].bin_prev = UINT32_MAX;

    if(allocator->bin_indices[bin_index] == UINT32_MAX)
    {
      allocator->used_bins[top_bin_index] &= ~(1 << leaf_bin_index);

      if(allocator->used_bins[top_bin_index] == 0)
        allocator->used_bins_top &= ~(1 << top_bin_index);
    }
  }

  allocator->free_nodes[allocator->free_offset++] = node_index;
  allocator->free_storage -= node->data_size;
}

/*
 * Create bookkeeping for a managed offset range of size bytes, initially all free.
 * One malloc holds the allocator, node records, and spare-node stack. The caller
 * owns the managed storage. max_allocations limits records for both used and free
 * blocks; one record describes the initial free block.
 * Returns NULL when max_allocations is zero or the metadata allocation fails.
 */
ftlsf_allocator *
ftlsf_create_allocator(_t(u32) size, _t(u32) max_allocations)
{
  if(max_allocations == 0)
    return NULL;

  _t(u64) nodes_offset = ftlsf_align_up(sizeof(ftlsf_allocator), sizeof(ftlsf_node));
  _t(u64) free_offset  = ftlsf_align_up(nodes_offset + max_allocations * sizeof(ftlsf_node), sizeof(_t(u32)));
  _t(u64) total_size   = free_offset + max_allocations * sizeof(_t(u32));

  ftlsf_allocator *allocator = malloc(total_size);
  if(allocator == NULL)
    return NULL;

  allocator->nodes           = (ftlsf_node *)((_t(u8) *)allocator + nodes_offset);
  allocator->free_nodes      = (_t(u32) *)((_t(u8) *)allocator + free_offset);
  allocator->max_allocations = max_allocations;
  allocator->size            = size;
  allocator->free_storage    = 0;
  allocator->used_bins_top   = 0;
  allocator->free_offset     = max_allocations;

  for(_t(u32) i = 0; i < FTLSF_NUM_TOP_BINS; ++i)
    allocator->used_bins[i] = 0;

  for(_t(u32) i = 0; i < FTLSF_NUM_LEAF_BINS; ++i)
    allocator->bin_indices[i] = UINT32_MAX;

  for(_t(u32) i = 0; i < max_allocations; ++i)
  {
    allocator->nodes[i].data_offset = 0;
    allocator->nodes[i].data_size   = 0;
    allocator->nodes[i].bin_next    = UINT32_MAX;
    allocator->nodes[i].bin_prev    = UINT32_MAX;
    allocator->nodes[i].bin_n_next  = UINT32_MAX;
    allocator->nodes[i].bin_n_prev  = UINT32_MAX;
    allocator->nodes[i].used        = 0;
  }

  for(_t(u32) i = 0; i < max_allocations; ++i)
    allocator->free_nodes[i] = max_allocations - i - 1;

  ftlsf_insert_node_into_bin(allocator, size, 0);

  return allocator;
}

/*
 * Release the allocator's bookkeeping and invalidate all handles belonging to it.
 * The managed storage remains owned by the caller. Passing NULL is allowed.
 */
void
ftlsf_destroy_allocator(ftlsf_allocator *allocator)
{
  free(allocator);
}

/*
 * Reserve size bytes from a suitable free bin and return their offset and node index.
 * Split any remainder into a free block, updating bins, neighbors, and free_storage.
 * A split requires a spare node; an exact fit reuses the existing record.
 * Returns UINT32_MAX in both fields if the size exceeds supported limits or a
 * suitable block or required spare node is unavailable. No offset alignment is added.
 * Requires a live allocator. A zero-size request can still consume a node record.
 */
ftlsf_allocation
ftlsf_allocate(ftlsf_allocator *allocator, _t(u64) size)
{
  if(size > UINT32_MAX || size > allocator->size)
  {
    return (ftlsf_allocation){.offset = UINT32_MAX, .meta = UINT32_MAX};
  }

  _t(u32) min_bin_index      = ftlsf_u32_to_float_round_up((_t(u32))size);
  _t(u32) min_top_bin_index  = min_bin_index >> FTLSF_TOP_BINS_INDEX_SHIFT;
  _t(u32) min_leaf_bin_index = min_bin_index & FTLSF_LEAF_BINS_INDEX_MASK;

  _t(u32) top_bin_index  = min_top_bin_index;
  _t(u32) leaf_bin_index = UINT32_MAX;

  if(allocator->used_bins_top & (1 << top_bin_index))
    leaf_bin_index = ftlsf_find_lowest_set_bit_after(allocator->used_bins[top_bin_index], min_leaf_bin_index);

  if(leaf_bin_index == UINT32_MAX)
  {
    top_bin_index = ftlsf_find_lowest_set_bit_after(allocator->used_bins_top, min_top_bin_index + 1);

    if(top_bin_index == UINT32_MAX)
      return (ftlsf_allocation){.offset = UINT32_MAX, .meta = UINT32_MAX};

    leaf_bin_index = ftlsf_tzcnt_nonzero(allocator->used_bins[top_bin_index]);
  }

  _t(u32) bin_index       = (top_bin_index << FTLSF_TOP_BINS_INDEX_SHIFT) | leaf_bin_index;
  _t(u32) node_index      = allocator->bin_indices[bin_index];
  ftlsf_node *node        = &allocator->nodes[node_index];
  _t(u32) node_total_size = node->data_size;
  _t(u32) remaining_size  = node_total_size - size;

  if(remaining_size != 0 && allocator->free_offset == 0)
    return (ftlsf_allocation){.offset = UINT32_MAX, .meta = UINT32_MAX};

  node->data_size                   = size;
  node->used                        = 1;
  allocator->bin_indices[bin_index] = node->bin_next;

  if(node->bin_next != UINT32_MAX)
    allocator->nodes[node->bin_next].bin_prev = UINT32_MAX;

  allocator->free_storage -= node_total_size;

  if(allocator->bin_indices[bin_index] == UINT32_MAX)
  {
    allocator->used_bins[top_bin_index] &= ~(1 << leaf_bin_index);

    if(allocator->used_bins[top_bin_index] == 0)
      allocator->used_bins_top &= ~(1 << top_bin_index);
  }

  if(remaining_size != 0)
  {
    _t(u32) new_node_index = ftlsf_insert_node_into_bin(allocator, remaining_size, node->data_offset + size);

    if(node->bin_n_next != UINT32_MAX)
      allocator->nodes[node->bin_n_next].bin_n_prev = new_node_index;

    allocator->nodes[new_node_index].bin_n_prev = node_index;
    allocator->nodes[new_node_index].bin_n_next = node->bin_n_next;
    node->bin_n_next                            = new_node_index;
  }

  return (ftlsf_allocation){.offset = node->data_offset, .meta = node_index};
}

/*
 * Release a live allocation using its metadata index and merge it with any free
 * physical neighbors on either side. Recycle the absorbed records, insert the
 * combined free block into its size bin, and reconnect the surrounding blocks.
 * Updates free_storage and the spare-node count. Requires a live allocator and
 * a successful allocation from that allocator that has not already been freed.
 */
void
ftlsf_free(ftlsf_allocator *allocator, ftlsf_allocation allocation)
{
  FTLSF_ASSERT(allocation.meta != UINT32_MAX);
  if(!allocator->nodes)
    return;

  _t(u32) node_index = allocation.meta;
  ftlsf_node *node   = &allocator->nodes[node_index];

  _t(u32) offset = node->data_offset;
  _t(u32) size   = node->data_size;

  if((node->bin_n_prev != UINT32_MAX) && (allocator->nodes[node->bin_n_prev].used == 0))
  {
    ftlsf_node *prev_node = &allocator->nodes[node->bin_n_prev];
    offset                = prev_node->data_offset;
    size += prev_node->data_size;

    ftlsf_remove_node_from_bin(allocator, node->bin_n_prev);

    FTLSF_ASSERT(prev_node->bin_n_next == node_index);
    node->bin_n_prev = prev_node->bin_n_prev;
  }

  if((node->bin_n_next != UINT32_MAX) && (allocator->nodes[node->bin_n_next].used == 0))
  {
    ftlsf_node *next_node = &allocator->nodes[node->bin_n_next];
    size += next_node->data_size;

    ftlsf_remove_node_from_bin(allocator, node->bin_n_next);

    FTLSF_ASSERT(next_node->bin_n_prev == node_index);
    node->bin_n_next = next_node->bin_n_next;
  }

  _t(u32) bin_n_next = node->bin_n_next;
  _t(u32) bin_n_prev = node->bin_n_prev;

  allocator->free_nodes[allocator->free_offset++] = node_index;

  _t(u32) combined_node_index = ftlsf_insert_node_into_bin(allocator, size, offset);

  if(bin_n_next != UINT32_MAX)
  {
    allocator->nodes[combined_node_index].bin_n_next = bin_n_next;
    allocator->nodes[bin_n_next].bin_n_prev          = combined_node_index;
  }
  if(bin_n_prev != UINT32_MAX)
  {
    allocator->nodes[combined_node_index].bin_n_prev = bin_n_prev;
    allocator->nodes[bin_n_prev].bin_n_next          = combined_node_index;
  }
}

#undef _t
#endif
