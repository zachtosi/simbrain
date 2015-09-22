#include <jni.h> 
#include <math.h>
#include "fastonebigheader.h"
#include "org_simbrain_special_projects_plasticity_proj_IPIFRule.h"

const float INV_SQRT_2PI = (float) 1.0f / sqrt(6.283185f);
const float NOISE = 1.0f;
const float N2 = NOISE * NOISE * 2.0f;
const float _3SIGMA = NOISE * 3.0f;
const float dz = 2 * _3SIGMA / 10.0f;

static inline float fxx1(float x) {
	if (x <= 0) return 0;
	return (x * 0.01f) / ((x * 0.01f) + 1.0f);
}

static inline float gConv(float x_0) 
{
	float sum = 0.0f;
	float z = -_3SIGMA;
	while (z <= _3SIGMA) {
		sum += (INV_SQRT_2PI/NOISE) * fasterexp(-(z * z)/N2) 
		* fxx1(x_0 - z) * dz;
		z += dz;
	}
	return sum;
}

JNIEXPORT jfloat JNICALL Java_org_simbrain_special_1projects_plasticity_1proj_IPIFRule_gaussConv(JNIEnv *env, jobject obj, jfloat x)
{
	return gConv((float) x);
}