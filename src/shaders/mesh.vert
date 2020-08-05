#version 450

// Vertex attributes
layout(location = 0) in vec3 inPosition;
layout(location = 1) in vec3 inNormal;
layout(location = 2) in vec4 inColor;
layout(location = 3) in uint inMaterialIndex;

layout(location = 0) out vec4 outColor;
layout(location = 1) out flat uint outMaterialIndex;

layout(binding = 0) uniform UniformBufferObject {
    mat4 view;
    mat4 proj;
} ubo;

float rand(vec2 uv)
{
    return fract(sin(dot(uv.xy ,vec2(12.9898, 78.233))) * 43758.5453);
}

void main()
{
	vec2 seed = vec2(inNormal.x * 1429.0 + inNormal.z * 1239.0, inNormal.y * 583.0 + inNormal.z * 3029.0);
	float randomFromNormal = rand(seed);
	vec4 extraColor = vec4(randomFromNormal * 0.05);

	outColor = inColor + extraColor;
    outMaterialIndex = inMaterialIndex;

    gl_Position = ubo.proj * ubo.view * vec4(inPosition, 1.0);
}
