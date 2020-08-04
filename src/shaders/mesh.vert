#version 450

// Vertex attributes
layout(location = 0) in vec3 inPosition;
layout(location = 1) in vec3 inNormal;
layout(location = 2) in uint inMaterialIndex;

layout(location = 0) out vec3 outColor;
layout(location = 1) out flat uint outMaterialIndex;

layout(binding = 0) uniform UniformBufferObject {
    mat4 view;
    mat4 proj;
} ubo;

void main()
{
	outColor = inNormal;
    outMaterialIndex = inMaterialIndex;

    gl_Position = ubo.proj * ubo.view * vec4(inPosition, 1.0);
}
