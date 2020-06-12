#include "GeometryHelper.h"

#include "Geometry.h"

#include "../Utilities/Maths.h"
#include "../Utilities/Vectors.h"

Geometry* GeometryHelper::CreateTriangle()
{
    std::vector<Vertex> vertices;

    vertices.push_back(Vertex(Point3<Metre>(Metre(0), Metre(0.5f), Metre(0)), Vector3<Real>(), Vector2<Real>()));
    vertices.push_back(Vertex(Point3<Metre>(Metre(0.5f), Metre(-0.5f), Metre(0)), Vector3<Real>(), Vector2<Real>()));
    vertices.push_back(Vertex(Point3<Metre>(Metre(-0.5f), Metre(-0.5f), Metre(0)), Vector3<Real>(), Vector2<Real>()));

    std::vector<uint32> indices;
    indices.push_back(0);
    indices.push_back(1);
    indices.push_back(2);
    Geometry* geom = Geometry::CreateGeometry("Triangle", std::move(vertices), std::move(indices));
    geom->setColor(Color::Red());
    return geom;
}


Geometry* GeometryHelper::CreateGrid(Metre width, Metre depth, uint32 m, uint32 n, float uRepeatRatio, float vRepeatRatio)
{
	std::vector<Vertex> vertices;
	std::vector<uint32> indices;

	// TP2 : À compléter
    // Vous n'utiliserez pas les textures pour ce TP, alors vous n'avez pas à considérer les paramètres float uRepeatRatio et float vRepeatRatio.
    // Un vertex contient sa position, son vecteur normal et ses coordonnées de texture. Pour ce TP, les coordonnées de texture seront (0,0).
    // Vous pouvez vous inspirer du cylindre pour faire la grille.
	
	for (int x = 0; x < m; x++) {
		for (int z = 0; z < n; z++) {
			Point3<Metre> p1 = Point3<Metre> (Metre(-width / 2 + x * (width/ m)), Metre(), Metre(-depth / 2 + z * (depth / n)));
			Point3<Metre> p2 = Point3<Metre>(Metre(-width / 2 + x * (width / m)), Metre(), Metre(-depth / 2 + (z + 1) * (depth / n)));
			Point3<Metre> p3 = Point3<Metre>(Metre(-width / 2 + (x + 1) * (width / m)), Metre(), Metre(-depth / 2 + (z + 1) * (depth / n)));
			Point3<Metre> p4 = Point3<Metre>(Metre(-width / 2 + (x + 1) * (width / m)), Metre(), Metre(-depth / 2 + z * (depth / n)));

			// Face 1
			Vector3<Real> normalF1V1 = (p3 - p1).normalized().crossProduct((p2 - p1).normalized());
			Vertex f1v1 = Vertex(p1, normalF1V1, Vector2<Real>());
			Vertex f1v2 = Vertex(p2, normalF1V1, Vector2<Real>());
			Vertex f1v3 = Vertex(p3, normalF1V1, Vector2<Real>());

			vertices.push_back(f1v1);
			vertices.push_back(f1v2);
			vertices.push_back(f1v3);

			indices.push_back(x * n * 6 + z * 6);
			indices.push_back((x * n * 6 + z * 6) + 2);
			indices.push_back((x * n * 6 + z * 6) + 1);

			Vector3<Real> normalF1V2 = (p4 - p1).normalized().crossProduct((p3 - p1).normalized());
			Vertex f1v4 = Vertex(p1, normalF1V2, Vector2<Real>());
			Vertex f1v5 = Vertex(p4, normalF1V2, Vector2<Real>());
			Vertex f1v6 = Vertex(p3, normalF1V2, Vector2<Real>());

			vertices.push_back(f1v4);
			vertices.push_back(f1v5);
			vertices.push_back(f1v6);

			indices.push_back((x * n * 6 + z * 6) + 3);
			indices.push_back((x * n * 6 + z * 6) + 4);
			indices.push_back((x * n * 6 + z * 6) + 5);
		}
	}
    
    Geometry* geom = Geometry::CreateGeometry("Grid", std::move(vertices), std::move(indices));
    return geom;
}

Geometry* GeometryHelper::CreateGrid(Metre width, Metre depth, uint32 m, uint32 n, const Color& color)
{
    Geometry* grid = GeometryHelper::CreateGrid(width, depth, m, n, 1.0f, 1.0f);
    grid->setColor(color);
    return grid;
}

Geometry* GeometryHelper::CreateGrid(Metre width, Metre depth, uint32 m, uint32 n, const ColorRGB& color)
{
    return CreateGrid(width, depth, m, n, Color(color));
}

Geometry* GeometryHelper::CreateBox(Metre width, Metre height, Metre depth, Real textureRatioX, Real textureRatioY)
{
	std::vector<Vertex> vertices;
	std::vector<uint32> indices;

	// Faite
	// TP2 : à compléter
	// Vous pouvez vous inspirer du tétrahèdre pour faire le cube. Le tétrahèdre remplace temporairement la sphère que vous aurez à faire au tp3.

	// Cube
	Point3<Metre> p1 = Point3<Metre>(Metre(-width/2), Metre(-height/2), Metre(-depth/2));
	Point3<Metre> p2 = Point3<Metre>(Metre(-width / 2), Metre(height/2), Metre(-depth / 2));
	Point3<Metre> p3 = Point3<Metre>(Metre(width/2), Metre(height/2), Metre(-depth / 2));
	Point3<Metre> p4 = Point3<Metre>(Metre(width/2), Metre(-height / 2), Metre(-depth / 2));
	Point3<Metre> p5 = Point3<Metre>(Metre(-width / 2), Metre(-height / 2), Metre(depth/2));
	Point3<Metre> p6 = Point3<Metre>(Metre(-width / 2), Metre(height/2), Metre(depth/2));
	Point3<Metre> p7 = Point3<Metre>(Metre(width/2), Metre(height/2), Metre(depth/2));
	Point3<Metre> p8 = Point3<Metre>(Metre(width/2), Metre(-height / 2), Metre(depth/2));

	// Face 1
	Vector3<Real> normalF1V1 = (p3 - p1).normalized().crossProduct((p2 - p1).normalized());
	Vertex f1v1 = Vertex(p1, normalF1V1, Vector2<Real>());
	Vertex f1v2 = Vertex(p2, normalF1V1, Vector2<Real>());
	Vertex f1v3 = Vertex(p3, normalF1V1, Vector2<Real>());

	vertices.push_back(f1v1);
	vertices.push_back(f1v2);
	vertices.push_back(f1v3);

	indices.push_back(0);
	indices.push_back(2);
	indices.push_back(1);

	Vector3<Real> normalF1V2 = (p4 - p1).normalized().crossProduct((p3 - p1).normalized());
	Vertex f1v4 = Vertex(p1, normalF1V2, Vector2<Real>());
	Vertex f1v5 = Vertex(p4, normalF1V2, Vector2<Real>());
	Vertex f1v6 = Vertex(p3, normalF1V2, Vector2<Real>());

	vertices.push_back(f1v4);
	vertices.push_back(f1v5);
	vertices.push_back(f1v6);

	indices.push_back(3);
	indices.push_back(4);
	indices.push_back(5);

	// Face 2
	Vector3<Real> normalF2V1 = (p4 - p1).normalized().crossProduct((p5 - p1).normalized());
	Vertex f2v1 = Vertex(p5, normalF2V1, Vector2<Real>());
	Vertex f2v2 = Vertex(p1, normalF2V1, Vector2<Real>());
	Vertex f2v3 = Vertex(p4, normalF2V1, Vector2<Real>());

	vertices.push_back(f2v1);
	vertices.push_back(f2v2);
	vertices.push_back(f2v3);

	indices.push_back(6);
	indices.push_back(7);
	indices.push_back(8);

	Vector3<Real> normalF2V2 = (p8 - p5).normalized().crossProduct((p4 - p5).normalized());
	Vertex f2v4 = Vertex(p5, normalF2V2, Vector2<Real>());
	Vertex f2v5 = Vertex(p8, normalF2V2, Vector2<Real>());
	Vertex f2v6 = Vertex(p4, normalF2V2, Vector2<Real>());

	vertices.push_back(f2v4);
	vertices.push_back(f2v5);
	vertices.push_back(f2v6);

	indices.push_back(9);
	indices.push_back(10);
	indices.push_back(11);

	// Face 3
	Vector3<Real> normalF3V1 = (p6 - p5).normalized().crossProduct((p7 - p5).normalized());
	Vertex f3v1 = Vertex(p5, normalF3V1, Vector2<Real>());
	Vertex f3v2 = Vertex(p6, normalF3V1, Vector2<Real>());
	Vertex f3v3 = Vertex(p7, normalF3V1, Vector2<Real>());

	vertices.push_back(f3v1);
	vertices.push_back(f3v2);
	vertices.push_back(f3v3);

	indices.push_back(12);
	indices.push_back(13);
	indices.push_back(14);
	
	Vector3<Real> normalF3V2 = (p8 - p5).normalized().crossProduct((p7 - p5).normalized());
	Vertex f3v4 = Vertex(p5, normalF3V2, Vector2<Real>());
	Vertex f3v5 = Vertex(p8, normalF3V2, Vector2<Real>());
	Vertex f3v6 = Vertex(p7, normalF3V2, Vector2<Real>());

	vertices.push_back(f3v4);
	vertices.push_back(f3v5);
	vertices.push_back(f3v6);

	indices.push_back(15);
	indices.push_back(16);
	indices.push_back(17);

	// Face 4
	Vector3<Real> normalF4V1 = (p7 - p6).normalized().crossProduct((p2 - p6).normalized());
	Vertex f4v1 = Vertex(p2, normalF4V1, Vector2<Real>());
	Vertex f4v2 = Vertex(p6, normalF4V1, Vector2<Real>());
	Vertex f4v3 = Vertex(p7, normalF4V1, Vector2<Real>());

	vertices.push_back(f4v1);
	vertices.push_back(f4v2);
	vertices.push_back(f4v3);

	indices.push_back(18);
	indices.push_back(19);
	indices.push_back(20);
	
	Vector3<Real> normalF4V2 = (p3 - p2).normalized().crossProduct((p7 - p2).normalized());
	Vertex f4v4 = Vertex(p2, normalF4V2, Vector2<Real>());
	Vertex f4v5 = Vertex(p3, normalF4V2, Vector2<Real>());
	Vertex f4v6 = Vertex(p7, normalF4V2, Vector2<Real>());

	vertices.push_back(f4v4);
	vertices.push_back(f4v5);
	vertices.push_back(f4v6);

	indices.push_back(21);
	indices.push_back(22);
	indices.push_back(23);

	// Face 5
	Vector3<Real> normalF5V1 = (p7 - p3).normalized().crossProduct((p4 - p3).normalized());
	Vertex f5v1 = Vertex(p4, normalF5V1, Vector2<Real>());
	Vertex f5v2 = Vertex(p3, normalF5V1, Vector2<Real>());
	Vertex f5v3 = Vertex(p7, normalF5V1, Vector2<Real>());

	vertices.push_back(f5v1);
	vertices.push_back(f5v2);
	vertices.push_back(f5v3);

	indices.push_back(24);
	indices.push_back(25);
	indices.push_back(26);

	
	Vector3<Real> normalF5V2 = (p8 - p4).normalized().crossProduct((p7 - p4).normalized());
	Vertex f5v4 = Vertex(p4, normalF5V2, Vector2<Real>());
	Vertex f5v5 = Vertex(p8, normalF5V2, Vector2<Real>());
	Vertex f5v6 = Vertex(p7, normalF5V2, Vector2<Real>());

	vertices.push_back(f5v4);
	vertices.push_back(f5v5);
	vertices.push_back(f5v6);

	indices.push_back(27);
	indices.push_back(28);
	indices.push_back(29);

	// Face 6
	Vector3<Real> normalF6V1 = (p2 - p6).normalized().crossProduct((p5 - p6).normalized());
	Vertex f6v1 = Vertex(p5, normalF6V1, Vector2<Real>());
	Vertex f6v2 = Vertex(p6, normalF6V1, Vector2<Real>());
	Vertex f6v3 = Vertex(p2, normalF6V1, Vector2<Real>());

	vertices.push_back(f6v1);
	vertices.push_back(f6v2);
	vertices.push_back(f6v3);

	indices.push_back(30);
	indices.push_back(31);
	indices.push_back(32);

	
	Vector3<Real> normalF6V2 = (p1 - p5).normalized().crossProduct((p2 - p5).normalized());
	Vertex f6v4 = Vertex(p5, normalF6V2, Vector2<Real>());
	Vertex f6v5 = Vertex(p1, normalF6V2, Vector2<Real>());
	Vertex f6v6 = Vertex(p2, normalF6V2, Vector2<Real>());

	vertices.push_back(f6v4);
	vertices.push_back(f6v5);
	vertices.push_back(f6v6);

	indices.push_back(33);
	indices.push_back(34);
	indices.push_back(35);


    Geometry* geom = Geometry::CreateGeometry("Cube", std::move(vertices), std::move(indices));
    return geom;
}

Geometry* GeometryHelper::CreateBox(Metre size, const Color& color)
{
    return CreateBox(size, size, size, color);
}

Geometry* GeometryHelper::CreateBox(Metre size, const ColorRGB& color)
{
    return CreateBox(size, Color(color));
}

Geometry* GeometryHelper::CreateBox(Metre width, Metre height, Metre depth)
{
    return CreateBox(width, height, depth, 1.0f, 1.0f);
}

Geometry* GeometryHelper::CreateBox(Metre width, Metre height, Metre depth, const Color& color)
{
    Geometry* box = CreateBox(width, height, depth, 1.0f, 1.0f);
    box->setColor(color);
    return box;
}

Geometry* GeometryHelper::CreateBox(Metre width, Metre height, Metre depth, const ColorRGB& color)
{
    return CreateBox(width, height, depth, Color(color));
}

Geometry* GeometryHelper::CreateSphere(Metre radius, uint32 slices, uint32 stacks)
{
    std::vector<Vertex> vertices;
    std::vector<uint32> indices;

	// TP3 : À compléter
	// Pour le TP3, vous devez remplacer le code du tétrahèdre suivant par celui d'une sphère.

	// Thetrahedron
	Point3<Metre> p1 = Point3<Metre>(-radius, Metre(-radius), Metre());
	Point3<Metre> p2 = Point3<Metre>(radius, Metre(-radius), Metre());
	Point3<Metre> p3 = Point3<Metre>(Metre(), Metre(radius), Metre(-radius));
	Point3<Metre> p4 = Point3<Metre>(Metre(), Metre(radius), Metre(radius));

	// Face 1
	Vector3<Real> normalF1 = (p4 - p1).normalized().crossProduct((p3 - p1).normalized());
	Vertex f1v1 = Vertex(p1, normalF1, Vector2<Real>());
	Vertex f1v2 = Vertex(p3, normalF1, Vector2<Real>());
	Vertex f1v3 = Vertex(p4, normalF1, Vector2<Real>());

	vertices.push_back(f1v1);
	vertices.push_back(f1v2);
	vertices.push_back(f1v3);

	indices.push_back(0);
	indices.push_back(2);
	indices.push_back(1);

	// Face 2
	Vector3<Real> normalF2 = (p3 - p2).normalized().crossProduct((p4 - p2).normalized());
	Vertex f2v1 = Vertex(p2, normalF2, Vector2<Real>());
	Vertex f2v2 = Vertex(p3, normalF2, Vector2<Real>());
	Vertex f2v3 = Vertex(p4, normalF2, Vector2<Real>());

	vertices.push_back(f2v1);
	vertices.push_back(f2v2);
	vertices.push_back(f2v3);

	indices.push_back(3);
	indices.push_back(4);
	indices.push_back(5);

	// Face 3
	Vector3<Real> normalF3 = (p1 - p4).normalized().crossProduct((p2 - p4).normalized());
	Vertex f3v1 = Vertex(p4, normalF3, Vector2<Real>());
	Vertex f3v2 = Vertex(p1, normalF3, Vector2<Real>());
	Vertex f3v3 = Vertex(p2, normalF3, Vector2<Real>());

	vertices.push_back(f3v1);
	vertices.push_back(f3v2);
	vertices.push_back(f3v3);

	indices.push_back(6);
	indices.push_back(7);
	indices.push_back(8);

	// Face 4
	Vector3<Real> normalF4 = (p2 - p3).normalized().crossProduct((p1 - p3).normalized());
	Vertex f4v1 = Vertex(p3, normalF4, Vector2<Real>());
	Vertex f4v2 = Vertex(p1, normalF4, Vector2<Real>());
	Vertex f4v3 = Vertex(p2, normalF4, Vector2<Real>());

	vertices.push_back(f4v1);
	vertices.push_back(f4v2);
	vertices.push_back(f4v3);

	indices.push_back(9);
	indices.push_back(11);
	indices.push_back(10);
	// Fin du code à compléter du TP3

    Geometry* geom = Geometry::CreateGeometry("Sphere", std::move(vertices), std::move(indices));
    return geom;
}

Geometry* GeometryHelper::CreateSphere(Metre radius, uint32 slices, uint32 stacks, const Color& color)
{
    Geometry* sphere = CreateSphere(radius, slices, stacks);
    sphere->setColor(color);
    return sphere;
}

Geometry* GeometryHelper::CreateSphere(Metre radius, uint32 slices, uint32 stacks, const ColorRGB& color)
{
    return CreateSphere(radius, slices, stacks, Color(color));
}
            
Geometry* GeometryHelper::CreateCylinder(Metre topRadius, Metre bottomRadius, Metre height, uint32 sliceCount, uint32 stackCount)
{
    std::vector<Vertex> vertices;
    std::vector<uint32> indices;
    
	// Build Stacks.
    Metre stackHeight = height / stackCount;

    // Amount to increment radius as we move up each stack level from bottom to top.
    Metre radiusStep = (topRadius - bottomRadius) / stackCount;

    uint32 ringCount = stackCount + 1;

    // Compute vertices for each stack ring starting at the bottom and moving up.
    for (uint32 i = 0; i < ringCount; ++i)
    {
        Metre y = -0.5f * height + i * stackHeight;
        Metre r = bottomRadius + i * radiusStep;

        // vertices of ring
        Degree dTheta = 2.0f * Degree(180) / sliceCount;
        for (uint32 j = 0; j <= sliceCount; ++j)
        {
            float c = Maths::Cos(j * dTheta);
            float s = Maths::Sin(j * dTheta);

            Point3<Metre> position = Point3<Metre>(r * c, y, r * s);

            // Cylinder can be parameterized as follows, where we introduce v
            // parameter that goes in the same direction as the v tex-coord
            // so that the bitangent goes in the same direction as the v tex-coord.
            //   Let r0 be the bottom radius and let r1 be the top radius.
            //   y(v) = h - hv for v in [0,1].
            //   r(v) = r1 + (r0-r1)v
            //
            //   x(t, v) = r(v)*cos(t)
            //   y(t, v) = h - hv
            //   z(t, v) = r(v)*sin(t)
            // 
            //  dx/dt = -r(v)*sin(t)
            //  dy/dt = 0
            //  dz/dt = +r(v)*cos(t)
            //
            //  dx/dv = (r0-r1)*cos(t)
            //  dy/dv = -h
            //  dz/dv = (r0-r1)*sin(t)

            // This is unit length.
            Vector3<Real> tangentU = Vector3<Real>(-s, Real(), c);

            Metre dr = bottomRadius - topRadius;
            Vector3<Metre> bitangent = Vector3<Metre>(dr * c, -height, dr * s);

            Vector3<Real> normal = tangentU.crossProduct(bitangent).normalized();

            Vertex v = Vertex(position, normal, Vector2<Real>(Real(float(j)) / sliceCount, 1.0f - Real(float(i)) / stackCount));

            vertices.push_back(v);
        }
    }

    // Add one because we duplicate the first and last vertex per ring
    // since the texture coordinates are different.
    uint32 ringVertexCount = sliceCount + 1;

    // Compute indices for each stack.
    for (uint32 i = 0; i < stackCount; ++i)
    {
        for (uint32 j = 0; j < sliceCount; ++j)
        {
            indices.push_back(i * ringVertexCount + j);
            indices.push_back((i + 1) * ringVertexCount + j);
			indices.push_back((i + 1) * ringVertexCount + j + 1);

            indices.push_back(i * ringVertexCount + j);
            indices.push_back((i + 1) * ringVertexCount + j + 1);
			indices.push_back(i * ringVertexCount + j + 1);
        }
    }

    // Build cylinder top cap
    uint32 baseIndexTop = (uint32)vertices.size();

    Metre yTop = 0.5f * height;
    Degree dThetaTop = 2.0f * Degree(180) / sliceCount;

    // Duplicate cap ring vertices because the texture coordinates and normals differ.
    for (uint32 i = 0; i <= sliceCount; ++i)
    {
        Metre x = topRadius * Maths::Cos(i * dThetaTop);
        Metre z = topRadius * Maths::Sin(i * dThetaTop);

        // Scale down by the height to try and make top cap texture coord area
        // proportional to base.
        float u = x / height + 0.5f;
        float v = z / height + 0.5f;

        Vertex vert = Vertex(Point3<Metre>(x, yTop, z), Vector3<Real>(Real(), Real(1), Real()), Vector2<Real>(u, v));
        vertices.push_back(vert);
    }

    // Cap center vertex.
    Vertex capVert = Vertex(Point3<Metre>(Metre(), yTop, Metre()), Vector3<Real>(Real(), Real(1), Real()), Vector2<Real>(0.5f, 0.5f));
    vertices.push_back(capVert);

    // Index of center vertex.
    uint32 centerIndexTop = (uint32)(vertices.size() - 1);

    for (uint32 i = 0; i < sliceCount; ++i)
    {
        indices.push_back(centerIndexTop);
        indices.push_back(baseIndexTop + i + 1);
		indices.push_back(baseIndexTop + i);
    }

    //Build cylinder bottom cap
    uint32 baseIndexBottom = (uint32)vertices.size();
    Metre yBottom = -0.5f * height;

    // vertices of ring
    Degree dThetaBottom = Degree(360) / sliceCount;
    for (uint32 i = 0; i <= sliceCount; ++i)
    {
        Metre x = bottomRadius * Maths::Cos(i * dThetaBottom);
        Metre z = bottomRadius * Maths::Sin(i * dThetaBottom);

        // Scale down by the height to try and make top cap texture coord area
        // proportional to base.
        float u = x / height + 0.5f;
        float v = z / height + 0.5f;

        Vertex vert = Vertex(Point3<Metre>(x, yBottom, z), Vector3<Real>(Real(), Real(-1), Real()), Vector2<Real>(u, v));
        vertices.push_back(vert);
    }

    // Cap center vertex.
    capVert = Vertex(Point3<Metre>(Metre(), yBottom, Metre()), Vector3<Real>(Real(), Real(-1), Real()), Vector2<Real>(0.5f, 0.5f));
    vertices.push_back(capVert);

    // Cache the index of center vertex.
    uint32 centerIndexBottom = (uint32)(vertices.size() - 1);

    for (uint32 i = 0; i < sliceCount; ++i)
    {
        indices.push_back(centerIndexBottom);
        indices.push_back(baseIndexBottom + i);
		indices.push_back(baseIndexBottom + i + 1);
    }

    Geometry* geom = Geometry::CreateGeometry("Cylinder", std::move(vertices), std::move(indices));
    return geom;
}

Geometry* GeometryHelper::CreateCylinder(Metre topRadius, Metre bottomRadius, Metre height, uint32 sliceCount, uint32 stackCount, const Color& color)
{
    Geometry* cylinder = CreateCylinder(topRadius, bottomRadius, height, sliceCount, stackCount);
    cylinder->setColor(color);
    return cylinder;
}

Geometry* GeometryHelper::CreateCylinder(Metre topRadius, Metre bottomRadius, Metre height, uint32 sliceCount, uint32 stackCount, const ColorRGB& color)
{
    return CreateCylinder(topRadius, bottomRadius, height, sliceCount, stackCount, Color(color));
}

// Rings : Nombre d'anneau
// Sides : Nombre de points pour une anneau
// Radius : Distance du centre à l'anneau
// ringRadius : Rayon du cercle interne autour de l'anneau
Geometry* GeometryHelper::CreateTorus(Metre radius, Metre ringRadius, uint32 sides, uint32 rings, Color* color)
{
    std::vector<Vertex> vertices;
    std::vector<uint32> indices;

	// TP3 : Bonus à compléter
    // Remplacer le code du tétrahèdre suivant par celui d'un tore (anneau)
    // Thetrahedron
    Point3<Metre> p1 = Point3<Metre>(-radius, Metre(-radius), Metre());
    Point3<Metre> p2 = Point3<Metre>(radius, Metre(-radius), Metre());
    Point3<Metre> p3 = Point3<Metre>(Metre(), Metre(radius), Metre(-radius));
    Point3<Metre> p4 = Point3<Metre>(Metre(), Metre(radius), Metre(radius));

    // Face 1
    Vector3<Real> normalF1 = (p4 - p1).normalized().crossProduct((p3 - p1).normalized());
    Vertex f1v1 = Vertex(p1, normalF1, Vector2<Real>());
    Vertex f1v2 = Vertex(p3, normalF1, Vector2<Real>());
    Vertex f1v3 = Vertex(p4, normalF1, Vector2<Real>());

    vertices.push_back(f1v1);
    vertices.push_back(f1v2);
    vertices.push_back(f1v3);

    indices.push_back(0);
    indices.push_back(2);
    indices.push_back(1);

    // Face 2
    Vector3<Real> normalF2 = (p3 - p2).normalized().crossProduct((p4 - p2).normalized());
    Vertex f2v1 = Vertex(p2, normalF2, Vector2<Real>());
    Vertex f2v2 = Vertex(p3, normalF2, Vector2<Real>());
    Vertex f2v3 = Vertex(p4, normalF2, Vector2<Real>());

    vertices.push_back(f2v1);
    vertices.push_back(f2v2);
    vertices.push_back(f2v3);

    indices.push_back(3);
    indices.push_back(4);
    indices.push_back(5);

    // Face 3
    Vector3<Real> normalF3 = (p1 - p4).normalized().crossProduct((p2 - p4).normalized());
    Vertex f3v1 = Vertex(p4, normalF3, Vector2<Real>());
    Vertex f3v2 = Vertex(p1, normalF3, Vector2<Real>());
    Vertex f3v3 = Vertex(p2, normalF3, Vector2<Real>());

    vertices.push_back(f3v1);
    vertices.push_back(f3v2);
    vertices.push_back(f3v3);

    indices.push_back(6);
    indices.push_back(7);
    indices.push_back(8);

    // Face 4
    Vector3<Real> normalF4 = (p2 - p3).normalized().crossProduct((p1 - p3).normalized());
    Vertex f4v1 = Vertex(p3, normalF4, Vector2<Real>());
    Vertex f4v2 = Vertex(p1, normalF4, Vector2<Real>());
    Vertex f4v3 = Vertex(p2, normalF4, Vector2<Real>());

    vertices.push_back(f4v1);
    vertices.push_back(f4v2);
    vertices.push_back(f4v3);

    indices.push_back(9);
    indices.push_back(11);
    indices.push_back(10);
	// Fin du bonus...

    Geometry* geom = Geometry::CreateGeometry("Torus", std::move(vertices), std::move(indices));
    if (color != nullptr)
    {
        geom->setColor(*color);
    }
    return geom;
}

Geometry* GeometryHelper::CreateRevolutionSurface(const std::vector<Point2<Metre>>& slicePoint, uint32 precision)
{
	std::vector<Vertex> vertices;
	std::vector<uint32> indices;

	// TP3 : À compléter
	
	Geometry* geom = Geometry::CreateGeometry("Revolution", std::move(vertices), std::move(indices));
	geom->updateNormals();
	return geom;
}

            