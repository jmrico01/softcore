#version 450

layout(location = 0) in vec3 inColor;
layout(location = 1) in flat uint inMaterialIndex;

layout(location = 0) out vec4 outColor;
layout(location = 1) out uint outMaterialIndex;

void main()
{
	outColor = vec4(inColor, 1.0);
	outMaterialIndex = inMaterialIndex;
}
