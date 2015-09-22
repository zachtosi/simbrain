#include <jni.h> 
#include <math.h>
#include "fastonebigheader.h"

const float INV_SQRT_2PI = (float) 1.0f / sqrt(6.283185f);

static inline float fxx1(float x) {
	if (x <= 0) return 0;
	return (x * 0.01) / ((x * 0.01) + 1);
}

static inline float gConv(float x_0, float noise) 
{
	float sum = 0.0f;
	float z = -noise;
	float dz = noise / 10.0f;
	float n2 = 2 * noise * noise;
	while (z <= noise) {
		sum += (INV_SQRT_2PI/noise) * fasterexp(-(z * z)/n2) 
		* fxx1(x_0 - z) * dz;
		z += dz;
	}
	return sum;
}

JNIEXPORT jfloat JNICALL Java_org_simbrain_special_1projects_plasticity_1proj_ModifiedBCMRule_gaussConv (JNIEnv *env, jobject obj, jfloat x, jfloat noise)
{
	return gConv((float) x, (float) (3 * noise));
}