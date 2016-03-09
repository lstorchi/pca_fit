#ifndef _PCAFITTER_VERSION_INC
#define _PCAFITTER_VERSION_INC

/*! lstorchi: Versioning stuff 
 * 
 * The first version where this was added after rel_0_0_0. 
 *      this is usefull to track versioning expecially when 
 *      distributiong some test exe
 * */
#define PCAFITTER_MAJOR_VER 0
#define PCAFITTER_MINOR_VER 7
#define PCAFITTER_PATCH_VER 0

#define PCAFITTER_VERSION ( PCAFITTER_MAJOR_VER * 10000 \
        + PCAFITTER_MINOR_VER * 100 \
        + PCAFITTER_PATCH_VER )

#endif
