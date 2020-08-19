#version 450

layout(location = 0) in vec2 inUv;

layout(location = 0) out vec4 outColor;

layout(binding = 0) uniform sampler2D rasterizedColor;
layout(binding = 1) uniform sampler2D raytracedColor;

// TODO can use raytraced alpha channel for something

void main()
{
	vec2 pixelSize = 1.0 / vec2(textureSize(raytracedColor, 0));
	float nearbyBlend = 0.0;

	// vec4 colorRasterized = texture(rasterizedColor, inUv);
	// vec3 colorRasterizedAlphaMult = colorRasterized.rgb * colorRasterized.a;

	vec4 colorRasterized = texture(rasterizedColor, inUv);
	colorRasterized += texture(rasterizedColor, inUv + vec2( pixelSize.x, 0.0)) * nearbyBlend;
	colorRasterized += texture(rasterizedColor, inUv + vec2(0.0,  pixelSize.y)) * nearbyBlend;
	colorRasterized += texture(rasterizedColor, inUv + vec2(0.0, -pixelSize.y)) * nearbyBlend;
	colorRasterized += texture(rasterizedColor, inUv + vec2(-pixelSize.x, 0.0)) * nearbyBlend;

	vec4 colorRaytraced = texture(raytracedColor, inUv);
	colorRaytraced += texture(raytracedColor, inUv + vec2( pixelSize.x, 0.0)) * nearbyBlend;
	colorRaytraced += texture(raytracedColor, inUv + vec2(0.0,  pixelSize.y)) * nearbyBlend;
	colorRaytraced += texture(raytracedColor, inUv + vec2(0.0, -pixelSize.y)) * nearbyBlend;
	colorRaytraced += texture(raytracedColor, inUv + vec2(-pixelSize.x, 0.0)) * nearbyBlend;

	//vec3 colorBlended = colorRaytraced;
	//vec3 colorBlended = colorRasterized;
	float gpuWeight = colorRasterized.a;
	float cpuWeight = 1.0 - gpuWeight;
	vec3 colorBlended = colorRaytraced.rgb * cpuWeight + colorRasterized.rgb * gpuWeight;
	//vec3 colorBlended = mix(colorRaytraced, colorRasterized.rgb, colorRasterized.a);

    outColor = vec4(colorBlended, 1.0);
}
