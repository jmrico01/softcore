#version 450

layout(location = 0) in vec2 inUv;

layout(location = 0) out vec4 outColor;

layout(binding = 0) uniform sampler2D rasterizedColor;
layout(binding = 1) uniform sampler2D raytracedColor;

void main()
{
	// TODO can use alpha channels for something (esp raytraced)

	// TODO blend these better. alpha = 1 on rasterized should mean no raytraced component (lerp?)
	vec4 colorRasterized = texture(rasterizedColor, inUv);
	vec3 colorRasterizedAlphaMult = colorRasterized.rgb * colorRasterized.a;

	vec4 colorRaytraced = texture(raytracedColor, inUv);

    outColor = vec4(colorRasterizedAlphaMult + colorRaytraced.rgb, 1.0);
}
