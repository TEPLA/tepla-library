#ifndef RDTSC_H_
#define RDTSC_H_

unsigned long long rdtsc()
{
	unsigned long long ret;
#ifdef __x86_64__
	__asm__ volatile ("rdtsc;shlq $32,%%rdx;orq %%rdx,%%rax":"=a"(ret)::"%rdx");
#else
	__asm__ volatile ("rdtsc" : "=A" (ret));
#endif
	return ret;
}

#endif /* RDTSC_H_ */
