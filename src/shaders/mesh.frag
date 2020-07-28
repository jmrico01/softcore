#version 450

layout(location = 0) in vec3 inNormal;
layout(location = 1) in vec3 inColor;

layout(location = 0) out vec4 outColor;

void main()
{
	vec3 lightAmbientColor = vec3(0.05, 0.05, 0.05);

	float lightDirIntensity = 0.8f;
	vec3 lightDirColor = vec3(1.0, 1.0, 1.0);
	vec3 lightDir1 = normalize(vec3(0.2, 0.6, -1.0));
	vec3 lightDir2 = normalize(vec3(-0.3, -0.1, 0.2));

	vec3 negNormal = -inNormal;
	float dot1 = max(0.0, dot(negNormal, lightDir1));
	float dot2 = max(0.0, dot(negNormal, lightDir2));

	vec3 totalLightDirColor = lightDirIntensity * lightDirColor * dot1 + lightDirIntensity * lightDirColor * dot2;

	vec3 totalLightColor = lightAmbientColor + totalLightDirColor;
	totalLightColor = clamp(totalLightColor, vec3(0.0, 0.0, 0.0), vec3(1.0, 1.0, 1.0));

	float gamma = 2.2;
	vec3 colorGammaCorrected = pow(inColor * totalLightColor, vec3(gamma));

	outColor = vec4(colorGammaCorrected, 1.0);
}
