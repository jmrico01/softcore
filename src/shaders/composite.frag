#version 450

struct Material {
	vec3 albedo;
	vec3 emissionColor;
	float smoothness;
	float emission;
};

layout(location = 0) in vec2 inUv;

layout(location = 0) out vec4 outColor;

layout(binding = 0) uniform sampler2D rasterizedColor;
layout(binding = 1) uniform sampler2D raytracedColor;

layout(binding = 2) uniform sampler2D triangleGeometry;
layout(binding = 3) uniform isampler2D triangleMaterials;

/*
layout(binding = 4) uniform UniformBufferObject {
	Material materials[8];
	vec3 filmTopLeft;
	vec3 filmUnitOffsetX;
	vec3 filmUnitOffsetY;
	vec3 cameraPos;
	uint numTriangles;
} ubo;

vec3 RaycastColor(vec3 rayOrigin, vec3 rayDir)
{
	uint triangleGeometryAtlasSize = 2048;
	uint triangleMaterialsAtlasSize = triangleGeometryAtlasSize / 4;

	float EPSILON = 0.000001;
	uint bounces = 4;
	float minDist = 0.1;
	float maxDist = 20.0;

	float intensity = 1.0;
	vec3 color = vec3(0.0, 0.0, 0.0);
	for (uint i = 0; i < bounces; i++) {
		uint hitMaterialIndex = 8;
		vec3 hitNormal;
		float hitDist = maxDist;

		for (uint j = 0; j < ubo.numTriangles; j++) {
			uint ind = j * 4;
			ivec2 baseInd = ivec2(ind % triangleGeometryAtlasSize, ind / triangleGeometryAtlasSize);
			vec3 triangleNormal = texelFetch(triangleGeometry, baseInd, 0).xyz;
			float dotNegRayNormal = dot(rayDir, triangleNormal);
			if (dotNegRayNormal > 0.0) {
				//continue;
			}

			vec3 a = texelFetch(triangleGeometry, baseInd + ivec2(1, 0), 0).xyz;
			vec3 b = texelFetch(triangleGeometry, baseInd + ivec2(2, 0), 0).xyz;
			vec3 c = texelFetch(triangleGeometry, baseInd + ivec2(3, 0), 0).xyz;

			// ray-triangle intersection
			vec3 ab = b - a;
			vec3 ac = c - a;
			vec3 h = cross(rayDir, ac);
			float x = dot(ab, h);
			if (x > -EPSILON && x < EPSILON) {
				continue;
			}

			float f = 1.0 / x;
			vec3 s = rayOrigin - a;
			float u = f * dot(s, h);
			if (u < 0.0 || u > 1.0) {
				continue;
			}

			vec3 q = cross(s, ab);
			float v = f * dot(rayDir, q);
			if (v < 0.0f || u + v > 1.0) {
				continue;
			}

			//return vec3(1.0, 1.0, 1.0);
			float t = f * dot(ac, q);
			if (t > minDist && t < hitDist) {
				ivec2 materialTexelInd = ivec2(j % triangleMaterialsAtlasSize, j / triangleMaterialsAtlasSize);
				hitMaterialIndex = texelFetch(triangleMaterials, materialTexelInd, 0).x;
				hitNormal = triangleNormal;
				hitDist = t;
			}
		}

		if (hitMaterialIndex == 8) {
			break;
		}

		Material hitMaterial = ubo.materials[hitMaterialIndex];
		if (hitMaterial.emission > 0.0) {
			color = hitMaterial.emission * hitMaterial.emissionColor * intensity;
			break;
		}
		else {
			intensity *= 0.5f;
			rayOrigin += rayDir * hitDist;

			vec3 pureBounce = rayDir - 2.0f * dot(rayDir, hitNormal) * hitNormal;
			rayDir = pureBounce;
		}
	}

	return color;
}
*/

// TODO can use raytraced alpha channel for something

void main()
{
	vec2 pixelSize = 1.0 / vec2(textureSize(raytracedColor, 0));

	vec4 colorRasterized = texture(rasterizedColor, inUv);
	vec3 colorRasterizedAlphaMult = colorRasterized.rgb * colorRasterized.a;

	float nearbyBlend = 0.4;
	vec3 colorRaytraced = texture(raytracedColor, inUv).rgb;
	colorRaytraced += texture(raytracedColor, inUv + vec2( pixelSize.x, 0.0)).rgb * nearbyBlend;
	colorRaytraced += texture(raytracedColor, inUv + vec2(0.0,  pixelSize.y)).rgb * nearbyBlend;
	colorRaytraced += texture(raytracedColor, inUv + vec2(0.0, -pixelSize.y)).rgb * nearbyBlend;
	colorRaytraced += texture(raytracedColor, inUv + vec2(-pixelSize.x, 0.0)).rgb * nearbyBlend;

	vec3 colorBlended = colorRaytraced;// mix(colorRaytraced, colorRasterized.rgb, colorRasterized.a);
    outColor = vec4(colorBlended, 1.0);

/*
	vec3 filmOffsetX = ubo.filmUnitOffsetX * gl_FragCoord.x;
	vec3 filmOffsetY = ubo.filmUnitOffsetY * gl_FragCoord.y;
	vec3 filmPos = ubo.filmTopLeft + filmOffsetX + filmOffsetY;
	vec3 rayDir = normalize(filmPos - ubo.cameraPos);

	//outColor = vec4(RaycastColor(ubo.cameraPos, rayDir), 1.0);
*/
}
