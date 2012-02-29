/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2009-2011 by Synge Todo <wistaria@comp-phys.org>,
*                            Haruhiko Matsuo <halm@looper.t.u-tokyo.ac.jp>
*
* This software is published under the ALPS Application License; you
* can use, redistribute it and/or modify it under the terms of the
* license, either version 1 or (at your option) any later version.
* 
* You should have received a copy of the ALPS Application License
* along with this software; see the file LICENSE. If not, the license
* is also available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#ifndef LOOPER_ATOMIC_IMPL_H
#define LOOPER_ATOMIC_IMPL_H

//
// compare_and_swap
//

#if defined(__linux__) && defined(__x86_64__) && defined(__FCC_VERSION) && !defined(_GNU_SOURCE)

// compare_and_swap assembler code written by Haruhiko Matsuo <halm@looper.t.u-tokyo.ac.jp>

// argument 1: %rdi (int& variable)
// argument 2: %esi (int  oldval)
// argument 3: %edx (int  newval)

asm(".align  16");
asm(".section        .gnu.linkonce.t.compare_and_swap,\"ax\",@progbits");
asm(".weak   compare_and_swap");
asm("compare_and_swap:");
asm("  movl    %esi, %eax");      // copy oldval to %eax
asm("  lock");
asm("  cmpxchgl %edx, (%rdi)");   // compare oldval(%eax) and *variable
asm("  sete %dl");                // set %dl value of ZF
asm("  xorl %eax, %eax");         // set %eax 0 (x86 common practice)
asm("  testb %dl, %dl");          // %dl && %dl. If 0, ZF=1 elseif 1, ZF=0
asm("  setne %al");               // retrun value (%al is a part of %eax)
asm("  ret");
asm(".size   compare_and_swap,.-compare_and_swap");
asm(".type   compare_and_swap,@function");


#elif defined(__linux__) && defined(__sparc) && defined(__FCC_VERSION) && !defined(_GNU_SOURCE)

asm(".section \".text\"");
asm(".align 4");
asm(".global atomic_cas_uint");
asm(".type atomic_cas_uint, #function");
asm("atomic_cas_uint:");
asm("      cas   [%o0],%o1,%o2");
asm("      retl");
asm("      mov   %o2,%o0");
asm(".size  atomic_cas_uint, (.-atomic_cas_uint)");

#endif

#endif // LOOPER_ATOMIC_IMPL_H
