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

	vec3 colorRasterized = texture(rasterizedColor, inUv).rgb;
	colorRasterized += texture(rasterizedColor, inUv + vec2( pixelSize.x, 0.0)).rgb * nearbyBlend;
	colorRasterized += texture(rasterizedColor, inUv + vec2(0.0,  pixelSize.y)).rgb * nearbyBlend;
	colorRasterized += texture(rasterizedColor, inUv + vec2(0.0, -pixelSize.y)).rgb * nearbyBlend;
	colorRasterized += texture(rasterizedColor, inUv + vec2(-pixelSize.x, 0.0)).rgb * nearbyBlend;

	vec3 colorRaytraced = texture(raytracedColor, inUv).rgb;
	colorRaytraced += texture(raytracedColor, inUv + vec2( pixelSize.x, 0.0)).rgb * nearbyBlend;
	colorRaytraced += texture(raytracedColor, inUv + vec2(0.0,  pixelSize.y)).rgb * nearbyBlend;
	colorRaytraced += texture(raytracedColor, inUv + vec2(0.0, -pixelSize.y)).rgb * nearbyBlend;
	colorRaytraced += texture(raytracedColor, inUv + vec2(-pixelSize.x, 0.0)).rgb * nearbyBlend;

	vec3 colorBlended = colorRasterized;// mix(colorRaytraced, colorRasterized.rgb, colorRasterized.a);
    outColor = vec4(colorBlended, 1.0);
}
