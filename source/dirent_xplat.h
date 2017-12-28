#ifndef DIRENT_XPLAT
#define DIRENT_XPLAT

/**
 * Header for cross-platform support for dirent
 *
 * Based on discussion:
 * http://stackoverflow.com/questions/612097/how-can-i-get-the-list-of-files-in-a-directory-using-c-or-c
 */

#ifdef _WIN32
#   include "dirent_win.h"
#else
#   include <dirent.h>
#endif // _WIN32

#endif // DIRENT_XPLAT
