#include"zkps.hpp"
using namespace std;
uint256 ADD(uint256 a, uint256 b, uint256 q)
{
	uint256 temp = a + b;//a<2^381 & b<2^381
	if(temp>=q)
		temp=temp-q;
	return temp;
}

uint256 SUB(uint256 a, uint256 b, uint256 q)
{
	uint256 temp;
	if(a>=b){
		temp = a - b;
	} else {
		temp = q - ( b - a );
	}
	return temp;
}
uint384 ADD0(uint384 a, uint384 b, const uint384 q)
{
	uint384 temp = a + b;//a<2^381 & b<2^381
	if(temp>=q)
		temp=temp-q;
	return temp;
}
uint384 SUB0(uint384 a, uint384 b, const uint384 q)
{
	uint384 temp;
	if(a>=b){
		temp = a - b;
	} else {
		temp = q - ( b - a );
	}
	return temp;
}
uint256 montgomery_reduce(uint512 t, const uint256 m, const uint256 inv)
{
	uint256 t_low = t(255, 0);
	uint512 temp=karatsuba_256(t_low, inv);
//	uint512 temp=t_low*inv;
	uint256 temp_low=temp(255,0);
	uint512 x;
	x=t+karatsuba_256(temp_low,m); //overflow warning!
//	x=t+temp_low*m;
	uint256 res=x(511,256);
	if(res>=m)
	{
		res=res-m;
	}
	return res;
}

ap_uint<512> karatsuba_256(ap_uint<256> x, ap_uint<256> y)
{
	if (x == 0 || y == 0) return 0;
    ap_uint<127> x_low = x(127, 0);
    ap_uint<127> x_high = x(255, 128);
    ap_uint<127> y_low = y(127, 0);
    ap_uint<127> y_high = y(255, 128);

    ap_uint<129> z1_left = x_low + x_high;
    ap_uint<129> z1_right = y_low + y_high;
    ap_uint<128> low1 = z1_left(127,0);
    ap_uint<128> low2 = z1_right(127,0);
    ap_uint<512> z0 = karatsuba_128(x_low, y_low);
    ap_uint<512> z2 = karatsuba_128(x_high, y_high);
    ap_uint<256> z1t = karatsuba_128(low1, low2);
    ap_uint<512> temp=z1_left[128]*z1_right[128];
    ap_uint<512> z1 = z1t+(ap_uint<256>(z1_left[128]*low2)<<128)+(ap_uint<256>(z1_right[128]*low1)<<128)+(temp<<256);
    ap_uint<512> res = (z2 << 256) + ((z1-z0-z2) << 128) + z0;
    return res;
}
uint384 montgomery_reduce0(uint768 t, const uint384& m, const uint384& inv)
{
	uint384 t_low = t(383, 0);
	uint768 temp=karatsuba_384(t_low, inv);
//	uint768 temp=t_low*inv;
	uint384 temp_low=temp(383,0);
	uint768 x;
	x=t+karatsuba_384(temp_low,m); //overflow warning!
//	x=t+temp_low*m;
	uint384 res=x(767,384);
	if(res>=m)
	{
		res=res-m;
	}
	return res;
}
uint384 MUL0(uint384 a, uint384 b, const uint384 q, const uint384 INV)
{
  uint768 temp=karatsuba_384(a,b);
//  uint768 temp=a*b;
  uint384 res=montgomery_reduce0(temp,q,INV);
  return res;
}


ap_uint<768> karatsuba_384(ap_uint<384> x, ap_uint<384> y)
{
#pragma HLS PIPELINE
    if (x == 0 || y == 0) return 0;

    ap_uint<128> x0 = x(127,0);
    ap_uint<128> x1 = x(255,128);
    ap_uint<128> x2 = x(383,256);

    ap_uint<128> y0 = y(127,0);
    ap_uint<128> y1 = y(255,128);
    ap_uint<128> y2 = y(383,256);

    ap_uint<256> z0 = karatsuba_128(x0,y0);
    ap_uint<256> z1 = karatsuba_128(x1,y1);
    ap_uint<256> z2 = karatsuba_128(x2,y2);
    ap_uint<256> z3 = karatsuba_128(x0 + x1,y0 + y1);
    ap_uint<256> z4 = karatsuba_128(x1 + x2,y1 + y2);
    ap_uint<256> z5 = karatsuba_128(x0 + x2,y0 + y2);

    ap_uint<256> p1 = z3 - z0 - z1;          // a0��b1 + a1��b0
    ap_uint<256> p3 = z4 - z1 - z2;          // a1��b2 + a2��b1
    ap_uint<256> p2 = z5 - z0 - z2 + z1;     // a0��b2 + a1��b1 + a2��b0

    ap_uint<768> result = (ap_uint<768>)z0 + (ap_uint<768>(p1) << 128) + (ap_uint<768>(p2) << 256) + (ap_uint<768>(p3) << 384) + (ap_uint<768>(z2) << 512);
    return result;
}


ap_uint<256> karatsuba_128(ap_uint<128> x, ap_uint<128> y)
{
	if (x == 0 || y == 0) return 0;
    ap_uint<64> x_low = x(63, 0);
    ap_uint<64> x_high = x(127, 64);
    ap_uint<64> y_low = y(63, 0);
    ap_uint<64> y_high = y(127, 64);

    ap_uint<65> z1_left = x_low + x_high;
    ap_uint<65> z1_right = y_low + y_high;
    ap_uint<64> low1 = z1_left(63,0);
    ap_uint<64> low2 = z1_right(63,0);
    ap_uint<256> z0 = x_low*y_low;
    ap_uint<256> z2 = x_high*y_high;
    ap_uint<128> z1t = low1*low2;
//    ap_uint<256> z0 = karatsuba_64(x_low, y_low);
//    ap_uint<256> z2 = karatsuba_64(x_high, y_high);
//    ap_uint<128> z1t = karatsuba_64(low1, low2);
    ap_uint<256> temp=z1_left[64]*z1_right[64];
    ap_uint<256> z1 = z1t+(ap_uint<128>(z1_left[64]*low2)<<64)+(ap_uint<128>(z1_right[64]*low1)<<64)+(temp<<128);
    ap_uint<256> res = (z2 << 128) + ((z1-z0-z2) << 64) + z0;
    return res;
}

uint256 MUL(uint256 a,uint256 b, const uint256 m, const uint256 inv)
{
//	uint512 temp=a*b;
	uint512 temp=karatsuba_256(a,b);
	uint256 res=montgomery_reduce(temp,m,inv);
	return res;
}

void read_w(uint256* ws,data_stream* w_s1,data_stream* w_s2,data_stream* w_s3,data_stream* w_s4)
{
	read_w_label1:for(int i=0;i<R*LOGR/2;i++)
	{
#pragma HLS PIPELINE
		uint256 tempw=ws[i];
		w_s1->write(tempw);
		w_s2->write(tempw);
		w_s3->write(tempw);
		w_s4->write(tempw);
	}
}


void read_s(uint256* scalars,s_stream* s1 )
{
	for(uint32 i=0;i<S;i++)
	{
#pragma HLS PIPELINE II=57
		uint256 temps=scalars[i];
		inddata s_ary[WINDOW_NUM];
		for(int j=0;j<WINDOW_NUM-1;j++)
		{
#pragma HLS PIPELINE off
			ap_uint<WINDOW_SIZE> temp=temps((j+1)*WINDOW_SIZE-1,j*WINDOW_SIZE);
			if(temp[WINDOW_SIZE-1]==1)
			{
				temps+=ONE<<((j+1)*WINDOW_SIZE);
			}
			s_ary[j]=temp(WINDOW_SIZE-1,0);
		}
		s_ary[WINDOW_NUM-1]=temps(255,(WINDOW_NUM-1)*WINDOW_SIZE);
		for(int k=0;k<WINDOW_NUM;k++)
		{
			s1->write(s_ary[k]);
		}
	}
}

void read_coeffs(uint256* coeffs,data_stream* c_s,uint32 step)
{
	uint32 startadd=(step*R)>>1;
	read_coeffs_label1:for(int i=0;i<R/2;i++)
	{
#pragma HLS PIPELINE
		uint256 tempc=coeffs[startadd+i];
		c_s->write(tempc);
	}
}


void subcore(data_stream* s1,data_stream* s2,data_stream* s3)
{
	subcore_label0:for(int i=0;i<R;i++)
	{
#pragma HLS PIPELINE
		uint256 t1,t2,t3;
		t1=s1->read();
		t2=s2->read();
		t3=SUB(t1,t2,bls377_p);
		s3->write(t3);
	}
}
void mulcore(data_stream* s1,data_stream* s2,data_stream* s3)
{
	mulcore_label0:for(int i=0;i<R;i++)
	{
#pragma HLS PIPELINE
		uint256 t1,t2,t3;
		t1=s1->read();
		t2=s2->read();
		t3=MUL(t1,t2,bls377_p,bls377_i);
		s3->write(t3);
	}
}
void writeback(uint256* nums, data_stream* n_s,uint32 step)
{
	uint32 startadd=(step*R)>>1;
	for(int i=0;i<R/2;i++)
	{
#pragma HLS PIPELINE
		uint256 tep=n_s->read();
		nums[startadd+i]=tep;
	}
}

void lazy_read(uint256* p_1,uint256* p_2,uint256* p_3,uint256* p_4,uint256* p_5,uint256* p_6,
		data_stream0* x_s,data_stream0* y_s,data_stream0* z_s,data_stream0* t_s)
{
	for(int i=0;i<N;i++)
	{
#pragma HLS PIPELINE
		uint384 xt,yt,zt,tt;
		xt=(uint384(p_1[i])<<128)+p_5[i](127,0);
		yt=(uint384(p_2[i])<<128)+p_5[i](255,128);
		zt=(uint384(p_3[i])<<128)+p_6[i](127,0);
		tt=(uint384(p_4[i])<<128)+p_6[i](255,128);
		x_s->write(xt);
		y_s->write(yt);
		z_s->write(zt);
		t_s->write(tt);
	}
}

void lazy_write(uint256* p_1,uint256* p_2,uint256* p_3,uint256* p_4,uint256* p_5,uint256* p_6,data_stream0* x_s,data_stream0* y_s,data_stream0* z_s,data_stream0* t_s)
{
	for(int i=0;i<32;i++)
	{
#pragma HLS PIPELINE
		uint256 p1,p2,p3,p4,p5,p6;
		uint384 xt=x_s->read();
		uint384 yt=y_s->read();
		uint384 zt=z_s->read();
		uint384 tt=t_s->read();
		p1=xt(383,128);
		p2=yt(383,128);
		p3=zt(383,128);
		p4=tt(383,128);
		p5=uint256(xt(127,0))+uint256(yt(127,0)<<128);
		p6=uint256(zt(127,0))+uint256(tt(127,0)<<128);
		p_1[i]=p1;
		p_2[i]=p2;
		p_3[i]=p3;
		p_4[i]=p4;
		p_5[i]=p5;
		p_6[i]=p6;
	}
}

extern"C"{

void threemul(data_stream0* in1,data_stream0* in2,data_stream0* out)
{
	for(int i=0;i<(WINDOWS+N)*3;)
	{
#pragma HLS PIPELINE II=1
		if (!in1->empty() && !in2->empty()) {
			uint384 t1,t2,t3;
			t1=in1->read();
			t2=in2->read();
			t3=MUL0(t1,t2,bls377_q,bls377_i);
			out->write(t3);
			i++;
		}
	}
}

void padds1(data_stream0* x1_s,data_stream0* y1_s,data_stream0* t1_s,data_stream0* x2_s,data_stream0* y2_s,data_stream0* t2_s,data_stream0* out1,data_stream0* out2)
{
	int counter = 0;
	uint384 ymx1, ymx2, ypx1, ypx2;
	uint384 x1,y1,t1,x2,y2,t2;

	for(int i=0;i<(WINDOWS+N)*3;i++)
	{
#pragma HLS PIPELINE II=1
		if (counter == 0) {
			x1=x1_s->read();
			y1=y1_s->read();
			x2=x2_s->read();
			y2=y2_s->read();
			t1=t1_s->read();
			t2=t2_s->read();
			ymx1=SUB0(y1,x1,bls377_q);
			ymx2=SUB0(y2,x2,bls377_q);
			ypx1=ADD0(x1,y1,bls377_q);
			ypx2=ADD0(x2,y2,bls377_q);
		}

		if (counter == 0) {
			out1->write(ymx1);
			out2->write(ymx2);
		} else if (counter == 1) {
			out1->write(ypx1);
			out2->write(ypx2);
		} else {
			out1->write(t1);
			out2->write(t2);
		}
		counter = (counter == 2) ? 0 : counter + 1;
	}
}

void padds2(data_stream0* m1_s,data_stream0* z1_s,data_stream0* z2_s,data_stream0* e_s,data_stream0* h_s,data_stream0* out1,data_stream0* out2)
{
	int counter = 0;
	uint384 A, B, ttt, E, H, z1, z2;

	for(int i=0;i<(WINDOWS+N)*3;i++)
	{
#pragma HLS PIPELINE II=1
		if (counter == 0) {
			A=m1_s->read();
			z1=z1_s->read();
			z2=z2_s->read();
		} else if (counter == 1) {
			B=m1_s->read();
			E=SUB0(B,A,bls377_q);
			H=ADD0(B,A,bls377_q);
			e_s->write(E);
			h_s->write(H);
		} else {
			ttt=m1_s->read();
		}

		if (counter == 0) {
			out1->write(bls377_k);
			out2->write(ttt);
		} else if (counter == 1) {
			out1->write(z1);
			out2->write(z2);
		} else {
			out1->write(E);
			out2->write(H);
		}
		counter = (counter == 2) ? 0 : counter + 1;
	}
}

void padds3(data_stream0* m2_s,data_stream0* e_s,data_stream0* h_s,data_stream0* out1,data_stream0* out2,data_stream0* t_s)
{
	int counter = 0;
	uint384 C, ztz, rt, E, H, D, F, G;

	for(int i=0;i<(WINDOWS+N)*3;i++)
	{
#pragma HLS PIPELINE II=1
		if (counter == 0) {
			C=m2_s->read();
			E=e_s->read();
			H=h_s->read();
		} else if (counter == 1) {
			ztz=m2_s->read();
			D=ADD0(ztz,ztz,bls377_q);
			F=SUB0(D,C,bls377_q);
			G=ADD0(D,C,bls377_q);
		} else {
			rt=m2_s->read();
			t_s->write(rt);
		}

		if (counter == 0) {
			out1->write(E);
			out2->write(F);
		} else if (counter == 1) {
			out1->write(G);
			out2->write(H);
		} else {
			out1->write(F);
			out2->write(G);
		}
		counter = (counter == 2) ? 0 : counter + 1;
	}
}

void padds4(data_stream0* m3_s,data_stream0* x_s,data_stream0* y_s,data_stream0* z_s)
{
	int counter = 0;
	for(int i=0;i<(WINDOWS+N)*3;i++)
	{
#pragma HLS PIPELINE II=1
		uint384 x,y,z;
		if (counter == 0) {
			x=m3_s->read();
		} else if (counter == 1) {
			y=m3_s->read();
		} else {
			z=m3_s->read();
			x_s->write(x);
			y_s->write(y);
			z_s->write(z);
		}
		counter = (counter == 2) ? 0 : counter + 1;
	}
}

void selector1(data_stream0* xm1_s,data_stream0* ym1_s,data_stream0* zm1_s,data_stream0* tm1_s,data_stream0* xm2_s,data_stream0* ym2_s,data_stream0* zm2_s,data_stream0* tm2_s,
		data_stream0* xn1_s,data_stream0* yn1_s,data_stream0* zn1_s,data_stream0* tn1_s,data_stream0* xn2_s,data_stream0* yn2_s,data_stream0* zn2_s,data_stream0* tn2_s,
		data_stream0* xr1_s,data_stream0* yr1_s,data_stream0* zr1_s,data_stream0* tr1_s,data_stream0* xr2_s,data_stream0* yr2_s,data_stream0* zr2_s,data_stream0* tr2_s)
{
	epoint tempm1,tempm2,tempn1,tempn2;
	for(int i=0;i<WINDOWS+N;)
	{
#pragma HLS PIPELINE II=3
		if(!xm1_s->empty() && !ym1_s->empty() && !zm1_s->empty() && !tm1_s->empty() && !xm2_s->empty() && !ym2_s->empty() && !zm2_s->empty() && !tm2_s->empty())
		{
			tempm1.x=xm1_s->read();
			tempm1.y=ym1_s->read();
			tempm1.z=zm1_s->read();
			tempm1.t=tm1_s->read();
			tempm2.x=xm2_s->read();
			tempm2.y=ym2_s->read();
			tempm2.z=zm2_s->read();
			tempm2.t=tm2_s->read();
			xr1_s->write(tempm1.x);
			yr1_s->write(tempm1.y);
			zr1_s->write(tempm1.z);
			tr1_s->write(tempm1.t);
			xr2_s->write(tempm2.x);
			yr2_s->write(tempm2.y);
			zr2_s->write(tempm2.z);
			tr2_s->write(tempm2.t);
			i++;
		}
		else if(!xn1_s->empty() && !yn1_s->empty() && !zn1_s->empty() && !tn1_s->empty() && !xn2_s->empty() && !yn2_s->empty() && !zn2_s->empty() && !tn2_s->empty())
		{
			tempn1.x=xn1_s->read();
			tempn1.y=yn1_s->read();
			tempn1.z=zn1_s->read();
			tempn1.t=tn1_s->read();
			tempn2.x=xn2_s->read();
			tempn2.y=yn2_s->read();
			tempn2.z=zn2_s->read();
			tempn2.t=tn2_s->read();
			xr1_s->write(tempn1.x);
			yr1_s->write(tempn1.y);
			zr1_s->write(tempn1.z);
			tr1_s->write(tempn1.t);
			xr2_s->write(tempn2.x);
			yr2_s->write(tempn2.y);
			zr2_s->write(tempn2.z);
			tr2_s->write(tempn2.t);
			i++;
		}
	}
}

void selector2(data_stream0* xm1_s,data_stream0* ym1_s,data_stream0* zm1_s,data_stream0* tm1_s,data_stream0* xn1_s,data_stream0* yn1_s,data_stream0* zn1_s,data_stream0* tn1_s,
		data_stream0* xr1_s,data_stream0* yr1_s,data_stream0* zr1_s,data_stream0* tr1_s)
{
	epoint temp;
	for(int i=0;i<N;)
	{
#pragma HLS PIPELINE
		if(!xr1_s->empty() && !yr1_s->empty() && !zr1_s->empty() && !tr1_s->empty())
		{
			temp.x=xr1_s->read();
			temp.y=yr1_s->read();
			temp.z=zr1_s->read();
			temp.t=tr1_s->read();
			xm1_s->write(temp.x);
			ym1_s->write(temp.y);
			zm1_s->write(temp.z);
			tm1_s->write(temp.t);
			i++;
		}
	}
	for(int j=0;j<WINDOWS;)
	{
#pragma HLS PIPELINE
		if(!xr1_s->empty() && !yr1_s->empty() && !zr1_s->empty() && !tr1_s->empty())
		{
			temp.x=xr1_s->read();
			temp.y=yr1_s->read();
			temp.z=zr1_s->read();
			temp.t=tr1_s->read();
			xn1_s->write(temp.x);
			yn1_s->write(temp.y);
			zn1_s->write(temp.z);
			tn1_s->write(temp.t);
			j++;
		}
	}
}

void bucket_part(s_stream* i_s,data_stream0* x_s,data_stream0* y_s,data_stream0* z_s,data_stream0* t_s,
                data_stream0* x1_s,data_stream0* y1_s,data_stream0* z1_s,data_stream0* t1_s,data_stream0* x2_s,data_stream0* y2_s,data_stream0* z2_s,data_stream0* t2_s,
				data_stream0* xn1_s,data_stream0* yn1_s,data_stream0* zn1_s,data_stream0* tn1_s,data_stream0* xn2_s,data_stream0* yn2_s,data_stream0* zn2_s,data_stream0* tn2_s,
                data_stream0* xr_s,data_stream0* yr_s,data_stream0* zr_s,data_stream0* tr_s,data_stream0* xr2_s,data_stream0* yr2_s,data_stream0* zr2_s,data_stream0* tr2_s,
				data_stream0* x3_s,data_stream0* y3_s,data_stream0* z3_s,data_stream0* t3_s)
{
#pragma HLS INTERFACE mode=s_axilite port=return

#ifdef __SYNTHESIS__
	epoint buckets[WINDOWS/2];
#pragma HLS DEPENDENCE dependent=false type=inter variable=buckets
#else
	epoint* buckets = (epoint*)malloc(WINDOWS/2*sizeof(epoint));
#endif

	uint1 flags[WINDOWS/2];
#pragma HLS DEPENDENCE dependent=false type=inter variable=flags
	s_stream nor_s;
#pragma HLS STREAM depth=32 variable=nor_s
	inddata conf_s[64];
	epoint conf_p[64];
	int count=0;
	for (int i=0; i<WINDOWS/2; i++)
		flags[i] = 0;

	for(int i=0;i<N;i++)
	{
#pragma HLS PIPELINE II=3
		epoint ptemp;
		inddata index;
		if(!xr_s->empty() && !yr_s->empty() && !zr_s->empty() && !tr_s->empty())
		{
			epoint pt0;
			inddata id0;
			pt0.x=xr_s->read();
			pt0.y=yr_s->read();
			pt0.z=zr_s->read();
			pt0.t=tr_s->read();
			id0=nor_s.read();
			flags[id0]=0;
			buckets[id0]=pt0;
		}
		ptemp.x=x_s->read();
		ptemp.y=y_s->read();
		ptemp.z=z_s->read();
		ptemp.t=t_s->read();
		index=i_s->read();
		if(index==0)
			continue;
		if(index[WINDOW_SIZE-1]==1)
		{
			index=inddata((1<<WINDOW_SIZE)-1)-index;
			ptemp.y=bls377_q-ptemp.y;
		}
		if(flags[index]==0)
		{
			flags[index]=1;
			x1_s->write(buckets[index].x);
			y1_s->write(buckets[index].y);
			z1_s->write(buckets[index].z);
			t1_s->write(buckets[index].t);
			x2_s->write(ptemp.x);
			y2_s->write(ptemp.y);
			z2_s->write(ptemp.z);
			t2_s->write(ptemp.t);
			nor_s.write(index);
		}
		else
		{
			conf_p[count]=ptemp;
			conf_s[count]=index;
			count=count+1;
		}
	}

	for(int ii=0;ii<count;ii++)
	{
#pragma HLS PIPELINE II=3
		epoint ptemp=conf_p[ii];
		inddata index=conf_s[ii];
		nor_s.write(index);
		x1_s->write(buckets[index].x);
		y1_s->write(buckets[index].y);
		z1_s->write(buckets[index].z);
		t1_s->write(buckets[index].t);
		x2_s->write(ptemp.x);
		y2_s->write(ptemp.y);
		z2_s->write(ptemp.z);
		t2_s->write(ptemp.t);
	}

	while (!xr_s->empty()) {
	    epoint pt0;
	    inddata id0;
	    pt0.x = xr_s->read();
	    pt0.y = yr_s->read();
	    pt0.z = zr_s->read();
	    pt0.t = tr_s->read();
	    id0 = nor_s.read();
	    flags[id0] = 0;
	    buckets[id0] = pt0;
	}

	epoint tmp[32];
#pragma HLS DEPENDENCE dependent=false type=inter variable=tmp
	epoint tmp1[32];
#pragma HLS DEPENDENCE dependent=false type=inter variable=tmp1

	uint32 c1=0,c2=0;
    while(c1<WINDOWS || c2<WINDOWS)
    {
		if (!xr2_s->empty() && !yr2_s->empty() && !zr2_s->empty() && !tr2_s->empty()) {
			uint32 kget=c2;
			epoint eback;
			eback.x=xr2_s->read();
			eback.y=yr2_s->read();
			eback.z=zr2_s->read();
			eback.t=tr2_s->read();
			inddata cat=kget(4,0);
			if(kget[5]==0)
			{
				tmp[cat]=eback;
			}
			else
			{
				tmp1[cat]=eback;
			}
            c2++;
        }
        inddata mouse=c1(4,0);
		if(c1[5]==0)
		{
			uint32 rk=256*(mouse+1)-1;
			uint32 rabbit=rk-c1(13,6);
			xn1_s->write(buckets[rabbit].x);
			yn1_s->write(buckets[rabbit].y);
			zn1_s->write(buckets[rabbit].z);
			tn1_s->write(buckets[rabbit].t);
			xn2_s->write(tmp[mouse].x);
			yn2_s->write(tmp[mouse].y);
			zn2_s->write(tmp[mouse].z);
			tn2_s->write(tmp[mouse].t);
		}
		else
		{
			xn1_s->write(tmp1[mouse].x);
			yn1_s->write(tmp1[mouse].y);
			zn1_s->write(tmp1[mouse].z);
			tn1_s->write(tmp1[mouse].t);
			xn2_s->write(tmp[mouse].x);
			yn2_s->write(tmp[mouse].y);
			zn2_s->write(tmp[mouse].z);
			tn2_s->write(tmp[mouse].t);
		}
		c1++;
    }
	for(uint32 k=0;k<32;k++)
	{
		x3_s->write(tmp1[k].x);
		y3_s->write(tmp1[k].y);
		z3_s->write(tmp1[k].z);
		t3_s->write(tmp1[k].t);
	}
}

void ntt1(data_stream* left0,data_stream* right0,data_stream* left1,data_stream* left2, data_stream* right1,data_stream* right2,data_stream* w_s,data_stream* outa_s, data_stream* outb_s)
{
	for(uint32 j=0;j<R/2;j++)
	{
#pragma HLS PIPELINE
		uint256 ina,inb;
		ina=left0->read();
		inb=right0->read();
		uint256 w=w_s->read();
		uint256 outa=ADD(ina,inb,bls377_p);
		uint256 outb=SUB(ina,inb,bls377_p);
		outb=MUL(outb,w,bls377_p,bls377_i);
		outa_s->write(outa);
		outb_s->write(outb);
	}
	for(int i=1;i<LOGR;i++)
	{
#pragma HLS PIPELINE off
		for(uint32 j=0;j<R/2;j++)
		{
#pragma HLS PIPELINE
			uint256 ina,inb;
			if(!j[0])
			{
				ina=left1->read();
				inb=right1->read();
			}
			else
			{
				ina=left2->read();
				inb=right2->read();
			}
			uint256 w=w_s->read();
			uint256 outa=ADD(ina,inb,bls377_p);
			uint256 outb=SUB(ina,inb,bls377_p);
			outb=MUL(outb,w,bls377_p,bls377_i);
			outa_s->write(outa);
			outb_s->write(outb);
		}
	}
}
void ntt2(data_stream* outa_s, data_stream* outb_s, data_stream* left1,data_stream* left2,data_stream* right1,data_stream* right2,data_stream* c1,data_stream* c2)
{
	for(int i=0;i<LOGR-1;i++)
	{
#pragma HLS PIPELINE off
		for(uint32 j=0;j<R/2;j++)
		{
#pragma HLS PIPELINE
			uint256 ina=outa_s->read();
			uint256 inb=outb_s->read();
			if(j[i]==0)
			{
				left1->write(ina);
				left2->write(inb);
			}
			else
			{
				right1->write(ina);
				right2->write(inb);
			}
		}
	}
	for(uint32 k=0;k<R/2;k++)
	{
		uint256 a1=outa_s->read();
		uint256 b1=outb_s->read();
		c1->write(a1);
		c2->write(b1);
	}
}

void data_loader_dual(
    uint32 mode, uint32 step,
    uint256* mem0,  uint256* mem1,  uint256* mem2,  uint256* mem3,
    uint256* mem4,  uint256* mem5,  uint256* mem6,  uint256* mem7,
    uint256* mem8,  uint256* mem9,  uint256* mem10, uint256* mem11,
    uint256* mem12, uint256* mem13, uint256* mem14, uint256* mem15,
    uint256* mem16,
    data_stream* w_s1, data_stream* w_s2, data_stream* w_s3, data_stream* w_s4,
    data_stream* left0,  data_stream* right0,
    data_stream* left1,  data_stream* right1,
    data_stream* left2,  data_stream* right2,
    data_stream* left3,  data_stream* right3,
    data_stream* ca0, data_stream* cb0,
    data_stream* ca1, data_stream* cb1,
    data_stream* ca2, data_stream* cb2,
    data_stream* ca3, data_stream* cb3,
    s_stream*    s1,
    data_stream0* x_s,  data_stream0* y_s,  data_stream0* z_s,  data_stream0* t_s,
    data_stream0* xw_s, data_stream0* yw_s, data_stream0* zw_s, data_stream0* tw_s
) {
#pragma HLS INTERFACE mode=s_axilite port=mode bundle=control
#pragma HLS INTERFACE mode=s_axilite port=step bundle=control
#pragma HLS INTERFACE mode=s_axilite port=return bundle=control

#pragma HLS INTERFACE mode=m_axi bundle=gmem1  port=mem0
#pragma HLS INTERFACE mode=m_axi bundle=gmem2  port=mem1
#pragma HLS INTERFACE mode=m_axi bundle=gmem3  port=mem2
#pragma HLS INTERFACE mode=m_axi bundle=gmem4  port=mem3
#pragma HLS INTERFACE mode=m_axi bundle=gmem5  port=mem4
#pragma HLS INTERFACE mode=m_axi bundle=gmem6  port=mem5
#pragma HLS INTERFACE mode=m_axi bundle=gmem7  port=mem6
#pragma HLS INTERFACE mode=m_axi bundle=gmem8  port=mem7
#pragma HLS INTERFACE mode=m_axi bundle=gmem9  port=mem8
#pragma HLS INTERFACE mode=m_axi bundle=gmem10 port=mem9
#pragma HLS INTERFACE mode=m_axi bundle=gmem11 port=mem10
#pragma HLS INTERFACE mode=m_axi bundle=gmem12 port=mem11
#pragma HLS INTERFACE mode=m_axi bundle=gmem13 port=mem12
#pragma HLS INTERFACE mode=m_axi bundle=gmem14 port=mem13
#pragma HLS INTERFACE mode=m_axi bundle=gmem15 port=mem14
#pragma HLS INTERFACE mode=m_axi bundle=gmem16 port=mem15
#pragma HLS INTERFACE mode=m_axi bundle=gmem17 port=mem16

#pragma HLS DATAFLOW

    if (mode) {
        read_w(mem0, w_s1, w_s2, w_s3, w_s4);

        read_coeffs(mem1, left0,  step);
        read_coeffs(mem2, right0, step);
        read_coeffs(mem3, left1,  step);
        read_coeffs(mem4, right1, step);
        read_coeffs(mem5, left2,  step);
        read_coeffs(mem6, right2, step);
        read_coeffs(mem7, left3,  step);
        read_coeffs(mem8, right3, step);

        writeback(mem9,  ca0, step);
        writeback(mem10, cb0, step);
        writeback(mem11, ca1, step);
        writeback(mem12, cb1, step);
        writeback(mem13, ca2, step);
        writeback(mem14, cb2, step);
        writeback(mem15, ca3, step);
        writeback(mem16, cb3, step);
    } else {
        read_s(mem0, s1);
        lazy_read(mem1, mem2, mem3, mem4, mem5, mem6, x_s,  y_s,  z_s,  t_s);
        lazy_write(mem9, mem10, mem11, mem12, mem13, mem14, xw_s, yw_s, zw_s, tw_s);
    }
}
}
void data_loader_ntt(uint256* ws,uint256* coeffs0,uint256* coeffs1,uint256* coeffs2,uint256* coeffs3,uint256* coeffs4,uint256* coeffs5,uint256* coeffs6,uint256* coeffs7,
		uint256* cout0,uint256* cout1,uint256* cout2,uint256* cout3,uint256* cout4,uint256* cout5,uint256* cout6,uint256* cout7,uint32 step,
		data_stream* w_s1,data_stream* w_s2,data_stream* w_s3,data_stream* w_s4,
		data_stream* left0,data_stream* right0,data_stream* left1,data_stream* right1,data_stream* left2,data_stream* right2,data_stream* left3,data_stream* right3,
		data_stream* ca0,data_stream* cb0,data_stream* ca1,data_stream* cb1,data_stream* ca2,data_stream* cb2,data_stream* ca3,data_stream* cb3)
{
#pragma HLS INTERFACE mode=s_axilite port=return
#pragma HLS INTERFACE mode=m_axi bundle=gmem17 port=cout7
#pragma HLS INTERFACE mode=m_axi bundle=gmem16 port=cout6
#pragma HLS INTERFACE mode=m_axi bundle=gmem15 port=cout5
#pragma HLS INTERFACE mode=m_axi bundle=gmem14 port=cout4
#pragma HLS INTERFACE mode=m_axi bundle=gmem13 port=cout3
#pragma HLS INTERFACE mode=m_axi bundle=gmem12 port=cout2
#pragma HLS INTERFACE mode=m_axi bundle=gmem11 port=cout1
#pragma HLS INTERFACE mode=m_axi bundle=gmem10 port=cout0
#pragma HLS INTERFACE mode=m_axi bundle=gmem9 port=coeffs7
#pragma HLS INTERFACE mode=m_axi bundle=gmem8 port=coeffs6
#pragma HLS INTERFACE mode=m_axi bundle=gmem7 port=coeffs5
#pragma HLS INTERFACE mode=m_axi bundle=gmem6 port=coeffs4
#pragma HLS INTERFACE mode=m_axi bundle=gmem5 port=coeffs3
#pragma HLS INTERFACE mode=m_axi bundle=gmem4 port=coeffs2
#pragma HLS INTERFACE mode=m_axi bundle=gmem3 port=coeffs1
#pragma HLS INTERFACE mode=m_axi bundle=gmem2 port=coeffs0
#pragma HLS INTERFACE mode=m_axi bundle=gmem1 port=ws
#pragma HLS DATAFLOW
	read_w(ws,w_s1,w_s2,w_s3,w_s4);
	read_coeffs(coeffs0,left0,step);
	read_coeffs(coeffs1,right0,step);
	read_coeffs(coeffs2,left1,step);
	read_coeffs(coeffs3,right1,step);
	read_coeffs(coeffs4,left2,step);
	read_coeffs(coeffs5,right2,step);
	read_coeffs(coeffs6,left3,step);
	read_coeffs(coeffs7,right3,step);
	writeback(cout0,ca0,step);
	writeback(cout1,cb0,step);
	writeback(cout2,ca1,step);
	writeback(cout3,cb1,step);
	writeback(cout4,ca2,step);
	writeback(cout5,cb2,step);
	writeback(cout6,ca3,step);
	writeback(cout7,cb3,step);
}
void data_loader_msm(uint256* scalars,s_stream* s1,
		uint256* p_1,uint256* p_2,uint256* p_3,uint256* p_4,uint256* p_5,uint256* p_6,
		data_stream0* x_s,data_stream0* y_s,data_stream0* z_s,data_stream0* t_s,
		uint256* pw_1,uint256* pw_2,uint256* pw_3,uint256* pw_4,uint256* pw_5,uint256* pw_6,
		data_stream0* xw_s,data_stream0* yw_s,data_stream0* zw_s,data_stream0* tw_s)
{
#pragma HLS INTERFACE mode=s_axilite port=return
#pragma HLS INTERFACE mode=m_axi bundle=gmem13 port=scalars
#pragma HLS INTERFACE mode=m_axi bundle=gmem12 port=pw_6
#pragma HLS INTERFACE mode=m_axi bundle=gmem11 port=pw_5
#pragma HLS INTERFACE mode=m_axi bundle=gmem10 port=pw_4
#pragma HLS INTERFACE mode=m_axi bundle=gmem9 port=pw_3
#pragma HLS INTERFACE mode=m_axi bundle=gmem8 port=pw_2
#pragma HLS INTERFACE mode=m_axi bundle=gmem7 port=pw_1
#pragma HLS INTERFACE mode=m_axi bundle=gmem6 port=p_6
#pragma HLS INTERFACE mode=m_axi bundle=gmem5 port=p_5
#pragma HLS INTERFACE mode=m_axi bundle=gmem4 port=p_4
#pragma HLS INTERFACE mode=m_axi bundle=gmem3 port=p_3
#pragma HLS INTERFACE mode=m_axi bundle=gmem2 port=p_2
#pragma HLS INTERFACE mode=m_axi bundle=gmem1 port=p_1
#pragma HLS DATAFLOW

	read_s(scalars,s1);
	lazy_read(p_1,p_2,p_3,p_4,p_5,p_6,x_s,y_s,z_s,t_s);
	lazy_write(pw_1,pw_2,pw_3,pw_4,pw_5,pw_6,xw_s,yw_s,zw_s,tw_s);
}
//void midntt(uint256* ws,uint256* coeffs1,uint256* coeffs2,uint256* coeffs3,uint256* coeffs4,uint256* coeffs5,uint256* coeffs6,uint256* coeffs7,uint256* coeffs8,
//		uint256* cout1,uint256* cout2,uint256* cout3,uint256* cout4,uint256* cout5,uint256* cout6,uint256* cout7,uint256* cout8,uint32 step)
//{
//#pragma HLS INTERFACE mode=s_axilite port=step bundle=control
//#pragma HLS INTERFACE mode=m_axi bundle=gmem17 port=cout8
//#pragma HLS INTERFACE mode=m_axi bundle=gmem16 port=cout7
//#pragma HLS INTERFACE mode=m_axi bundle=gmem15 port=cout6
//#pragma HLS INTERFACE mode=m_axi bundle=gmem14 port=cout5
//#pragma HLS INTERFACE mode=m_axi bundle=gmem13 port=cout4
//#pragma HLS INTERFACE mode=m_axi bundle=gmem12 port=cout3
//#pragma HLS INTERFACE mode=m_axi bundle=gmem11 port=cout2
//#pragma HLS INTERFACE mode=m_axi bundle=gmem10 port=cout1
//#pragma HLS INTERFACE mode=m_axi bundle=gmem9 port=coeffs8
//#pragma HLS INTERFACE mode=m_axi bundle=gmem8 port=coeffs7
//#pragma HLS INTERFACE mode=m_axi bundle=gmem7 port=coeffs6
//#pragma HLS INTERFACE mode=m_axi bundle=gmem6 port=coeffs5
//#pragma HLS INTERFACE mode=m_axi bundle=gmem5 port=coeffs4
//#pragma HLS INTERFACE mode=m_axi bundle=gmem4 port=coeffs3
//#pragma HLS INTERFACE mode=m_axi bundle=gmem3 port=coeffs2
//#pragma HLS INTERFACE mode=m_axi bundle=gmem2 port=coeffs1
//#pragma HLS INTERFACE mode=m_axi bundle=gmem1 port=ws
//#pragma HLS DATAFLOW
//	data_stream w_s1,w_s2,w_s3,w_s4;
//#pragma HLS STREAM depth=16 variable=w_s4
//#pragma HLS STREAM depth=16 variable=w_s3
//#pragma HLS STREAM depth=16 variable=w_s2
//#pragma HLS STREAM depth=16 variable=w_s1
//	read_w(ws,&w_s1,&w_s2,&w_s3,&w_s4);
//	ntt0(&w_s1,coeffs1,coeffs2,step,cout1,cout2);
//	ntt0(&w_s2,coeffs3,coeffs4,step,cout3,cout4);
//	ntt0(&w_s3,coeffs5,coeffs6,step,cout5,cout6);
//	ntt0(&w_s4,coeffs7,coeffs8,step,cout7,cout8);
//}
//void step2(uint256* num1,uint256* num2,uint256* num3,uint32 step)
//{
//#pragma HLS INTERFACE mode=s_axilite port=step bundle=control
//#pragma HLS INTERFACE mode=m_axi bundle=gmem14 port=num1
//#pragma HLS INTERFACE mode=m_axi bundle=gmem15 port=num2
//#pragma HLS INTERFACE mode=m_axi bundle=gmem16 port=num3
//#pragma HLS DATAFLOW
//	data_stream s1,s2,s3;
//#pragma HLS STREAM depth=16 variable=s3
//#pragma HLS STREAM depth=16 variable=s2
//#pragma HLS STREAM depth=16 variable=s1
//	readnum(num1,&s1,step);
//	readnum(num2,&s2,step);
//	mulcore(&s1,&s2,&s3);
//	writeback(num3,&s3,step);
//}
//void hstep(uint256* num1,uint256* num2,uint256* num3,uint256* num4, uint32 step)
//{
//#pragma HLS INTERFACE mode=s_axilite port=step bundle=control
//#pragma HLS INTERFACE mode=m_axi bundle=gmem14 port=num1
//#pragma HLS INTERFACE mode=m_axi bundle=gmem15 port=num2
//#pragma HLS INTERFACE mode=m_axi bundle=gmem16 port=num3
//#pragma HLS INTERFACE mode=m_axi bundle=gmem17 port=num4
//#pragma HLS DATAFLOW
//	data_stream s1,s2,s3,s4,s5;
//#pragma HLS STREAM depth=16 variable=s5
//#pragma HLS STREAM depth=16 variable=s4
//#pragma HLS STREAM depth=16 variable=s3
//#pragma HLS STREAM depth=16 variable=s2
//#pragma HLS STREAM depth=16 variable=s1
//	readnum(num1,&s1,step);
//	readnum(num2,&s2,step);
//	readnum(num4,&s4,step);
//	subcore(&s1,&s4,&s5);
//	mulcore(&s5,&s2,&s3);
//	writeback(num3,&s3,step);
//}
void unipadd(data_stream0* x1_s,data_stream0* y1_s,data_stream0* z1_s,data_stream0* t1_s,data_stream0* x2_s,data_stream0* y2_s,data_stream0* z2_s,data_stream0* t2_s,
		data_stream0* xr_s,data_stream0* yr_s,data_stream0* zr_s,data_stream0* tr_s)
{
#pragma HLS DATAFLOW
	data_stream0 ina1,inb1,out1,e_s,h_s,ina2,inb2,out2,ina3,inb3,out3;
#pragma HLS STREAM depth=8 variable=out3
#pragma HLS STREAM depth=8 variable=inb3
#pragma HLS STREAM depth=8 variable=ina3
#pragma HLS STREAM depth=8 variable=out2
#pragma HLS STREAM depth=8 variable=inb2
#pragma HLS STREAM depth=8 variable=ina2
#pragma HLS STREAM depth=8 variable=h_s
#pragma HLS STREAM depth=8 variable=e_s
#pragma HLS STREAM depth=8 variable=out1
#pragma HLS STREAM depth=8 variable=inb1
#pragma HLS STREAM depth=8 variable=ina1
	padds1(x1_s,y1_s,t1_s,x2_s,y2_s,t2_s,&ina1,&inb1);
	threemul(&ina1,&inb1,&out1);
	padds2(&out1,z1_s,z2_s,&e_s,&h_s,&ina2,&inb2);
	threemul(&ina2,&inb2,&out2);
	padds3(&out2,&e_s,&h_s,&ina3,&inb3,tr_s);
	threemul(&ina3,&inb3,&out3);
	padds4(&out3,xr_s,yr_s,zr_s);
}

void ntt(uint256* ws,uint256* coeffs0,uint256* coeffs1,uint256* coeffs2,uint256* coeffs3,uint256* coeffs4,uint256* coeffs5,uint256* coeffs6,uint256* coeffs7,
		uint256* cout0,uint256* cout1,uint256* cout2,uint256* cout3,uint256* cout4,uint256* cout5,uint256* cout6,uint256* cout7,uint32 step)
{
#pragma HLS DATAFLOW
	data_stream w_s1,w_s2,w_s3,w_s4,left0,right0,left1,right1,left2,right2,left3,right3;
	data_stream left01,right01,left11,right11,left21,right21,left31,right31;
	data_stream outa0,outb0,outa1,outb1,outa2,outb2,outa3,outb3;
	data_stream left02,right02,left12,right12,left22,right22,left32,right32;
	data_stream ca0,cb0,ca1,cb1,ca2,cb2,ca3,cb3;
	ntt1(&left0,&right0,&left01,&left02,&right01,&right02,&w_s1,&outa0,&outb0);
	ntt2(&outa0,&outb0,&left01,&left02,&right01,&right02,&ca0,&cb0);
	ntt1(&left1,&right1,&left11,&left12,&right11,&right12,&w_s2,&outa1,&outb1);
	ntt2(&outa1,&outb1,&left11,&left12,&right11,&right12,&ca1,&cb1);
	ntt1(&left2,&right2,&left21,&left22,&right21,&right22,&w_s3,&outa2,&outb2);
	ntt2(&outa2,&outb2,&left21,&left22,&right21,&right22,&ca2,&cb2);
	ntt1(&left3,&right3,&left31,&left32,&right31,&right32,&w_s4,&outa3,&outb3);
	ntt2(&outa3,&outb3,&left31,&left32,&right31,&right32,&ca3,&cb3);
	data_loader_ntt(ws,coeffs0,coeffs1,coeffs2,coeffs3,coeffs4,coeffs5,coeffs6,coeffs7,cout0,cout1,cout2,cout3,cout4,cout5,cout6,cout7,step,
			&w_s1,&w_s2,&w_s3,&w_s4,&left0,&right0,&left1,&right1,&left2,&right2,&left3,&right3,&ca0,&cb0,&ca1,&cb1,&ca2,&cb2,&ca3,&cb3);
}
void msm(uint256* p_1,uint256* p_2,uint256* p_3,uint256* p_4,uint256* p_5,uint256* p_6,uint256* scalars,
		uint256* r_1,uint256* r_2,uint256* r_3,uint256* r_4,uint256* r_5,uint256* r_6)
{
#pragma HLS INTERFACE mode=s_axilite port=return
#pragma HLS INTERFACE mode=m_axi bundle=gmem13 port=scalars depth=1024
#pragma HLS INTERFACE mode=m_axi bundle=gmem12 port=r_6 depth=32
#pragma HLS INTERFACE mode=m_axi bundle=gmem11 port=r_5 depth=32
#pragma HLS INTERFACE mode=m_axi bundle=gmem10 port=r_4 depth=32
#pragma HLS INTERFACE mode=m_axi bundle=gmem9 port=r_3 depth=32
#pragma HLS INTERFACE mode=m_axi bundle=gmem8 port=r_2 depth=32
#pragma HLS INTERFACE mode=m_axi bundle=gmem7 port=r_1 depth=32
#pragma HLS INTERFACE mode=m_axi bundle=gmem6 port=p_6 depth=19456
#pragma HLS INTERFACE mode=m_axi bundle=gmem5 port=p_5 depth=19456
#pragma HLS INTERFACE mode=m_axi bundle=gmem4 port=p_4 depth=19456
#pragma HLS INTERFACE mode=m_axi bundle=gmem3 port=p_3 depth=19456
#pragma HLS INTERFACE mode=m_axi bundle=gmem2 port=p_2 depth=19456
#pragma HLS INTERFACE mode=m_axi bundle=gmem1 port=p_1 depth=19456
#pragma HLS DATAFLOW

	data_stream0 x1_s,y1_s,z1_s,t1_s,xm1_s,ym1_s,zm1_s,tm1_s,xm2_s,ym2_s,zm2_s,tm2_s;
#pragma HLS STREAM depth=16 variable=tm2_s
#pragma HLS STREAM depth=16 variable=zm2_s
#pragma HLS STREAM depth=16 variable=ym2_s
#pragma HLS STREAM depth=16 variable=xm2_s
#pragma HLS STREAM depth=16 variable=tm1_s
#pragma HLS STREAM depth=16 variable=zm1_s
#pragma HLS STREAM depth=16 variable=ym1_s
#pragma HLS STREAM depth=16 variable=xm1_s
#pragma HLS STREAM depth=16 variable=t1_s
#pragma HLS STREAM depth=16 variable=z1_s
#pragma HLS STREAM depth=16 variable=y1_s
#pragma HLS STREAM depth=16 variable=x1_s
	data_stream0 xn1_s,yn1_s,zn1_s,tn1_s,xn2_s,yn2_s,zn2_s,tn2_s,xc1_s,yc1_s,zc1_s,tc1_s,xc2_s,yc2_s,zc2_s,tc2_s;
#pragma HLS STREAM depth=16 variable=tc2_s
#pragma HLS STREAM depth=16 variable=zc2_s
#pragma HLS STREAM depth=16 variable=yc2_s
#pragma HLS STREAM depth=16 variable=xc2_s
#pragma HLS STREAM depth=16 variable=tc1_s
#pragma HLS STREAM depth=16 variable=zc1_s
#pragma HLS STREAM depth=16 variable=yc1_s
#pragma HLS STREAM depth=16 variable=xc1_s
#pragma HLS STREAM depth=16 variable=tn2_s
#pragma HLS STREAM depth=16 variable=zn2_s
#pragma HLS STREAM depth=16 variable=yn2_s
#pragma HLS STREAM depth=16 variable=xn2_s
#pragma HLS STREAM depth=16 variable=tn1_s
#pragma HLS STREAM depth=16 variable=zn1_s
#pragma HLS STREAM depth=16 variable=yn1_s
#pragma HLS STREAM depth=16 variable=xn1_s
	data_stream0 xi_s,yi_s,zi_s,ti_s,xj_s,yj_s,zj_s,tj_s;
#pragma HLS STREAM depth=16 variable=tj_s
#pragma HLS STREAM depth=16 variable=zj_s
#pragma HLS STREAM depth=16 variable=yj_s
#pragma HLS STREAM depth=16 variable=xj_s
#pragma HLS STREAM depth=16 variable=ti_s
#pragma HLS STREAM depth=16 variable=zi_s
#pragma HLS STREAM depth=16 variable=yi_s
#pragma HLS STREAM depth=16 variable=xi_s
	data_stream0 xr1_s,yr1_s,zr1_s,tr1_s,xr2_s,yr2_s,zr2_s,tr2_s;
#pragma HLS STREAM depth=16 variable=tr2_s
#pragma HLS STREAM depth=16 variable=zr2_s
#pragma HLS STREAM depth=16 variable=yr2_s
#pragma HLS STREAM depth=16 variable=xr2_s
#pragma HLS STREAM depth=16 variable=tr1_s
#pragma HLS STREAM depth=16 variable=zr1_s
#pragma HLS STREAM depth=16 variable=yr1_s
#pragma HLS STREAM depth=16 variable=xr1_s
	data_stream0 xr_s,yr_s,zr_s,tr_s,xf_s,yf_s,zf_s,tf_s;
#pragma HLS STREAM depth=16 variable=tf_s
#pragma HLS STREAM depth=16 variable=zf_s
#pragma HLS STREAM depth=16 variable=yf_s
#pragma HLS STREAM depth=16 variable=xf_s
#pragma HLS STREAM depth=16 variable=tr_s
#pragma HLS STREAM depth=16 variable=zr_s
#pragma HLS STREAM depth=16 variable=yr_s
#pragma HLS STREAM depth=16 variable=xr_s
	s_stream s1;
#pragma HLS STREAM depth=16 variable=s1

#ifdef __SYNTHESIS__
	bucket_part(&s1,&x1_s,&y1_s,&z1_s,&t1_s,&xm1_s,&ym1_s,&zm1_s,&tm1_s,&xm2_s,&ym2_s,&zm2_s,&tm2_s,&xn1_s,&yn1_s,&zn1_s,&tn1_s,&xn2_s,&yn2_s,&zn2_s,&tn2_s,&xr1_s,&yr1_s,&zr1_s,&tr1_s,&xr2_s,&yr2_s,&zr2_s,&tr2_s,&xf_s,&yf_s,&zf_s,&tf_s);
	selector1(&xm1_s,&ym1_s,&zm1_s,&tm1_s,&xm2_s,&ym2_s,&zm2_s,&tm2_s,&xn1_s,&yn1_s,&zn1_s,&tn1_s,&xn2_s,&yn2_s,&zn2_s,&tn2_s,&xc1_s,&yc1_s,&zc1_s,&tc1_s,&xc2_s,&yc2_s,&zc2_s,&tc2_s);
	unipadd(&xc1_s,&yc1_s,&zc1_s,&tc1_s,&xc2_s,&yc2_s,&zc2_s,&tc2_s,&xr_s,&yr_s,&zr_s,&tr_s);
	selector2(&xr1_s,&yr1_s,&zr1_s,&tr1_s,&xr2_s,&yr2_s,&zr2_s,&tr2_s,&xr_s,&yr_s,&zr_s,&tr_s);
	data_loader_msm(scalars,&s1,p_1,p_2,p_3,p_4,p_5,p_6,&x1_s,&y1_s,&z1_s,&t1_s,r_1,r_2,r_3,r_4,r_5,r_6,&xf_s,&yf_s,&zf_s,&tf_s);
#endif
}









