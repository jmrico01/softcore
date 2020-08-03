#version 450

layout(location = 0) in flat uint inMaterialIndex;

layout(location = 0) out uint outMaterialIndex;

void main()
{
	outMaterialIndex = inMaterialIndex;
}
