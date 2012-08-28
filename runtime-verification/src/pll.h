#pragma once

#define pll_e				0
#define pll_eb				1
#define pll_in				2
#define pll_inb				3
#define pll_mout			4
#define pll_moutb			5
#define pll_osc				6
#define pll_oscb			7
#define pll_out				8
#define pll_outb			9
#define pll_xvco_c			10
#define pll_xvco_s			11
#define pll_xvco_s_clip		12
#define pll_xpd_clip1		13
#define pll_xpd_clip2		14
#define pll_xpd_n1			15
#define pll_time			16


void setInitialPLLState(double* state);
