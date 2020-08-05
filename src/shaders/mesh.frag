#version 450

layout(location = 0) in vec4 inColor;
layout(location = 1) in flat uint inMaterialIndex;

layout(location = 0) out vec4 outColor;
layout(location = 1) out uint outMaterialIndex;

void main()
{
	outColor = inColor;
	outMaterialIndex = inMaterialIndex;
}
