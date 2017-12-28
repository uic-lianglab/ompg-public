#ifndef VLA_H
#define VLA_H

/**
 * Utilities for declaring variable length arrays (VLAs)
 */

/**
 * Visual Studio does not support VLAs. Can use alloca but, due to stack size
 * and non-POD data type concerns, will just use std::vector instead
 */
#ifdef _MSC_VER

// Note: alternative implementation could use alloca:
// https://msdn.microsoft.com/en-us/library/wb1s57t5(v=vs.120).aspx
// but would likely not work for non-POD types (constructors won't be called)

#include <vector>

#define DECLARE_VLA_1D(type_, name_, num_elem_) \
    std::vector<type_> name_(num_elem_)

#define DECLARE_VLA_2D(type_, name_, num_elem_1_, num_elem_2_) \
    std::vector<std::vector<type_> > name_(num_elem_1_, std::vector<type_>(num_elem_2_))

#define DECLARE_VLA_3D(type_, name_, num_elem_1_, num_elem_2_, num_elem_3_) \
    std::vector<std::vector< std::vector< type_> > > name_(num_elem_1_, std::vector<std::vector<type_> >(num_elem_2_, std::vector<type_>(num_elem_3_)))

#define SET_VLA_1D(type_, name_, num_elem_, val_) \
    memset(&name_[0], val_, num_elem_ * sizeof(type_))

#define SET_VLA_2D(type_, name_, num_elem_1_, num_elem_2_, val_) \
    do { \
        for (size_t i__ = 0; i__ < num_elem_1_; ++i__) { \
            memset(&(name_[i__][0]), val_, num_elem_2_ * sizeof(type_)); \
        } \
    } while(0)

#define SET_VLA_3D(type_, name_, num_elem_1_, num_elem_2_, num_elem_3_, val_) \
    do { \
        for (size_t i__ = 0; i__ < num_elem_1_; ++i__) { \
            for (size_t j__ = 0; j__ < num_elem_2_; ++j__) { \
                memset(&(name_[i__][j__][0]), val_, num_elem_3_ * sizeof(type_)); \
            } \
        } \
    } while(0)

#else

/**
 * Assuming GCC style compiler which supports C99
 */

#define DECLARE_VLA_1D(type_, name_, num_elem_) \
    type_ name_[num_elem_]

#define DECLARE_VLA_2D(type_, name_, num_elem_1_, num_elem_2_) \
    type_ name_[(num_elem_1_)][(num_elem_2_)]

#define DECLARE_VLA_3D(type_, name_, num_elem_1_, num_elem_2_, num_elem_3_) \
    type_ name_[(num_elem_1_)][(num_elem_2_)][(num_elem_3_)]

#define SET_VLA_1D(type_, name_, num_elem_, val_) \
    memset(&name_[0], val_, num_elem_ * sizeof(type_))

#define SET_VLA_2D(type_, name_, num_elem_1_, num_elem_2_, val_) \
    memset(&(name_[0][0]), val_, num_elem_1_ * num_elem_2_ * sizeof(type_))

#define SET_VLA_3D(type_, name_, num_elem_1_, num_elem_2_, num_elem_3_, val_) \
    memset(&(name_[0][0][0]), val_, num_elem_1_ * num_elem_2_ * num_elem_3_ * sizeof(type_))

#endif // _MSC_VER

#endif // VLA_H
