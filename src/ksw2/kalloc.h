#ifndef _KALLOC_H_
#define _KALLOC_H_

#include <stdlib.h>

#define km_size(x) (*(((size_t*)(x))-1) * sizeof(size_t))

#ifdef __cplusplus
extern "C" {
#endif

void *kmalloc_(void *km, size_t size);
void *krealloc_(void *km, void *ptr, size_t size);
void *kcalloc_(void *km, size_t count, size_t size);
void kfree_(void *km, void *ptr);

void *km_init_(void);
void km_destroy_(void *km);

void km_stat_(const void *km); // TODO: return numbers instead of print to stderr

#ifdef __cplusplus
}
#endif

#endif
