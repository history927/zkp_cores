#include<ap_int.h>
#include<hls_stream.h>
#define WINDOW_SIZE 14
#define WINDOWS (1<<14)
#define WINDOW_NUM 19
#define S 1024
#define N 19456
#define R 1024
#define LOGR 10
//N=S*WINDOW_NUM
//#define S (1<<22)
//#if (S==65536)
//	#define WINDOWS (1<<14)
//	#define WINDOW_SIZE 14
//	#define WINDOW_NUM 19
//	#define NS 1245184					//N=S*WINDOW_NUM
//	#define MULTIME 3735552				//MULTIME=NS*3
//	#define N 65536
//	#define R 256
//	#define LOGR 8
//#elif (S==524288)
//	#define WINDOWS (1<<14)
//	#define WINDOW_SIZE 14
//	#define WINDOW_NUM 19
//	#define NS 9961472					//N=S*WINDOW_NUM
//	#define MULTIME 14942208			//MULTIME=NS*3
//	#define N 524288
//	#define R 1024
//	#define LOGR 10
//#elif (S==4194304)
//	#define WINDOWS (1<<14)
//	#define WINDOW_SIZE 14
//	#define WINDOW_NUM 19
//	#define NS 79691776					//N=S*WINDOW_NUM
//	#define MULTIME 239075328			//MULTIME=NS*3
//	#define N 4194304
//	#define R 2048
//	#define LOGR 11
//#endif

typedef ap_uint<256> uint256;
typedef ap_uint<1> uint1;
typedef ap_uint<WINDOW_SIZE> inddata;
typedef ap_uint<32> uint32;
typedef ap_uint<512> uint512;
typedef ap_uint<384> uint384;
typedef ap_uint<768> uint768;
const uint384 bls381_q("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787");
const uint384 bls381_i("31812345668571394200522959936182890292651532161468729907319195075155855129916747767697203739041175453056836858347517");//for montgomery
const uint384 bls377_q("258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458177");
const uint384 bls377_i0("29496869032736653100683447689715864783909294504611102227839586092113932643754048458779746488862010077557953021345791");//for montgomery
const uint384 bls377_k("244536567197351118976972678317271058193963773829754279159068307164067353570771581460084726682472071493849921806358");
const uint256 bn256_p("21888242871839275222246405745257275088548364400416034343698204186575808495617");
const uint256 bls377_p("8444461749428370424248824938781546531375899335154063827935233455917409239041");
const uint256 bls377_i("47752251086953357377073236701509605140872345086634869599321669320666611974143");
const uint256 ONE=1;
ap_uint<256> karatsuba_128(ap_uint<128> x, ap_uint<128> y);
ap_uint<768> karatsuba_384(ap_uint<384> x, ap_uint<384> y);
uint256 montgomery_reduce(uint512 t, const uint256 m, const uint256 inv);
ap_uint<512> karatsuba_256(ap_uint<256> x, ap_uint<256> y);
struct epoint {
    uint384 x=0;
    uint384 y=1;
    uint384 z=0;
    uint384 t=0;
};
typedef hls::stream<epoint> ep_stream;
typedef hls::stream<uint384> data_stream0;
typedef hls::stream<inddata> s_stream;
typedef hls::stream<uint256> data_stream;
