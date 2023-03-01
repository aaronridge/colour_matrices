#include<stdio.h>
#include<cmath>
#include<cstring>

struct Primary
{
	float x;
	float y;
};

struct Colour_space
{
	Primary r;
	Primary g;
	Primary b;
	Primary w;
};

struct Vec3
{
	float v[3];
};

struct M33
{
	Vec3 r[3];
};

const M33 identity =
{
	Vec3{1.f, 0.f, 0.f},
	Vec3{0.f, 1.f, 0.f},
	Vec3{0.f, 0.f, 1.f}
};

float dot
(
	const Vec3& a,
	const Vec3& b
)
{
	return a.v[0] * b.v[0] + a.v[1] * b.v[1] + a.v[2] * b.v[2];
}

Vec3 mul
(
	const M33& m,
	const Vec3& v
)
{
	return Vec3
	{
		dot(m.r[0], v),
		dot(m.r[1], v),
		dot(m.r[2], v)
	};
}

Vec3 mul
(
	const Vec3& v,
	float scale
)
{
	return Vec3
	{
		v.v[0] * scale,
		v.v[1] * scale,
		v.v[2] * scale
	};
}

M33 invert
(
	const M33& m_in
)
{
	M33 res = identity;
	M33 m = m_in;
	Vec3 temp_v;

	// Arrange pivots
	for (int p = 0; p < 3; ++p)
	{
		float max_v = 0.f;
		int max_row = p;

		for (int r = p; r < 3; ++r)
		{
			const float v = fabsf(m.r[r].v[p]);
			if (v > max_v)
			{
				max_row = r;
				max_v = v;
			}
		}

		if (max_row == p)
		{
			continue;
		}

		temp_v = m.r[p];
		m.r[p] = m.r[max_row];
		m.r[max_row] = temp_v;

		temp_v = res.r[p];
		res.r[p] = res.r[max_row];
		res.r[max_row] = temp_v;
	}

	// Normalize rows
	for (int p = 0; p < 3; ++p)
	{
		if (m.r[p].v[p] == 0.f)
		{
			// No pivot, error
			return identity;
		}

		const float scale = 1.f / m.r[p].v[p];

		for (int c = 0; c < 3; ++c)
		{
			m.r[p].v[c] *= scale;
			res.r[p].v[c] *= scale;
		}

		// Eliminate from other rows
		for (int r = 0; r < 3; ++r)
		{
			if (r == p)
			{
				continue;
			}

			const float v = m.r[r].v[p];
			
			{
				m.r[r].v[0] -= v * m.r[p].v[0];
				m.r[r].v[1] -= v * m.r[p].v[1];
				m.r[r].v[2] -= v * m.r[p].v[2];
			}

			{
				res.r[r].v[0] -= v * res.r[p].v[0];
				res.r[r].v[1] -= v * res.r[p].v[1];
				res.r[r].v[2] -= v * res.r[p].v[2];
			}
		}
	}

	return res;
}

Vec3 primary_to_xyz
(
	const Primary& primary
)
{
	return Vec3
	{
		primary.x,
		primary.y,
		1.f - primary.x - primary.y
	};
}

M33 colour_space_to_xyz_matrix
(
	const Colour_space& colour_space,
	bool verbose
)
{
	if (verbose)
	{
		printf("Constructing matrix from colour space to XYZ:\n\n");
		printf("Primaries:\n"
			"Red: x = %.8f y = %.8f\n"
			"Green: x = %.8f y = %.8f\n"
			"Blue: x = %.8f y = %.8f\n"
			"White: x = %.8f y = %.8f\n\n",
			colour_space.r.x, colour_space.r.y,
			colour_space.g.x, colour_space.g.y,
			colour_space.b.x, colour_space.b.y,
			colour_space.w.x, colour_space.w.y);
	}

	// Convert primaries to XYZ vectors
	const Vec3 red = primary_to_xyz(colour_space.r);
	const Vec3 green = primary_to_xyz(colour_space.g);
	const Vec3 blue = primary_to_xyz(colour_space.b);
	const Vec3 white = primary_to_xyz(colour_space.w);

	if (verbose)
	{
		printf("Convert from xy to XYZ\n"
			"xy coordinate is the projection of XYZ onto X + Y + Z = 1 plane\n"
			"X = x\n"
			"Y = y\n"
			"Z = 1 - x - y\n\n");
		printf("Converting primaries and white point to XYZ vectors:\n"
			"Red = (%.8f, %.8f, %.8f)\n"
			"Green = (%.8f, %.8f, %.8f)\n"
			"Blue = (%.8f, %.8f, %.8f)\n"
			"White = (%.8f, %.8f, %.8f)\n\n",
			red.v[0], red.v[1], red.v[2],
			green.v[0], green.v[1], green.v[2],
			blue.v[0], blue.v[1], blue.v[2],
			white.v[0], white.v[1], white.v[2]
		);
	}

	const M33 base_matrix =
	{
		Vec3{red.v[0], green.v[0], blue.v[0]},
		Vec3{red.v[1], green.v[1], blue.v[1]},
		Vec3{red.v[2], green.v[2], blue.v[2]}
	};

	if (verbose)
	{
		printf("The XYZ primary vectors must be independently scaled such that they sum to white\n"
			"This means White = M * S\n"
			"Where M = [ Xr, Xg, Xb ]\n"
			"          [ Yr, Yg, Yb ]\n"
			"          [ Zr, Zg, Zb ]\n"
			"= [ %.8f, %.8f, %.8f ]\n"
			"  [ %.8f, %.8f, %.8f ]\n"
			"  [ %.8f, %.8f, %.8f ]\n\n"
			"For the primaries provided\n\n",
			base_matrix.r[0].v[0], base_matrix.r[0].v[1], base_matrix.r[0].v[2],
			base_matrix.r[1].v[0], base_matrix.r[1].v[1], base_matrix.r[1].v[2],
			base_matrix.r[2].v[0], base_matrix.r[2].v[1], base_matrix.r[2].v[2]
		);
	}

	// Invert basis matrix
	const M33 inverse_matrix = invert(base_matrix);

	if (verbose)
	{
		printf("To calculate S, we must multiply W by the inverse of M:\n"
			"S = InvM * W\n\n"
			"The inverse of M\n"
			"= [ %.8f, %.8f, %.8f ]\n"
			"  [ %.8f, %.8f, %.8f ]\n"
			"  [ %.8f, %.8f, %.8f ]\n\n"
			"For the primaries provided\n\n",
			inverse_matrix.r[0].v[0], inverse_matrix.r[0].v[1], inverse_matrix.r[0].v[2],
			inverse_matrix.r[1].v[0], inverse_matrix.r[1].v[1], inverse_matrix.r[1].v[2],
			inverse_matrix.r[2].v[0], inverse_matrix.r[2].v[1], inverse_matrix.r[2].v[2]
		);
	}

	// Multiply inverse basis by white point to calculate primary scales
	Vec3 scales = mul(inverse_matrix, white);

	if (verbose)
	{
		printf ("This gives S as:\n"
			"(%.8f, %.8f, %.8f)\n",
			scales.v[0], scales.v[1], scales.v[2]
		);
		printf("The white point is specified such that the luminance, or Yw, should equal 1\n"
			"This means we must divide S by Yw so that when our primaries are summed together the results Y value equals 1\n"
			"(%.8f, %.8f, %.8f) / %.8f = ",
			scales.v[0], scales.v[1], scales.v[2], white.v[1]
		);
	}

	// Normalize by luminance = white.y
	scales.v[0] /= white.v[1];
	scales.v[1] /= white.v[1];
	scales.v[2] /= white.v[1];

	if (verbose)
	{
		printf("(%.8f, %.8f, %.8f)\n",
			scales.v[0], scales.v[1], scales.v[2]
		);
	}

	const Vec3 scaled_red = mul(red, scales.v[0]);
	const Vec3 scaled_green = mul(green, scales.v[1]);
	const Vec3 scaled_blue = mul(blue, scales.v[2]);

	if (verbose)
	{
		printf("Scaling our primaries by S gives:\n"
			"Red = (%.8f, %.8f, %.8f)\n"
			"Green = (%.8f, %.8f, %.8f)\n"
			"Blue = (%.8f, %.8f, %.8f)\n\n",
			scaled_red.v[0], scaled_red.v[1], scaled_red.v[2], 
			scaled_green.v[0], scaled_green.v[1], scaled_green.v[2], 
			scaled_blue.v[0], scaled_blue.v[1], scaled_blue.v[2]
		);
	}

	return M33
	{
		Vec3{scaled_red.v[0], scaled_green.v[0], scaled_blue.v[0]},
		Vec3{scaled_red.v[1], scaled_green.v[1], scaled_blue.v[1]},
		Vec3{scaled_red.v[2], scaled_green.v[2], scaled_blue.v[2]}
	};
}

int main
(
	int argc,
	char** argv
)
{
	Colour_space space = { };

	bool found[4] = { false, false, false, false };

	bool verbose = false;

	int i = 0;
	while (i < argc)
	{
		if (strcmp(argv[i], "--v") == 0)
		{
			verbose = true;
			++i;
			continue;
		}

		if (strcmp(argv[i], "--help") == 0)
		{
			printf("Colour space tool\n"
				"Enter a colour space as xy primaries\n"
				"Use option -r [x] [y] for red primary\n"
				"Use option -g [x] [y] for green primary\n"
				"Use option -b [x] [y] for blue primary\n"
				"Use option -w [x] [y] for white primary\n"
			);

			return 0;
		}

		if (strcmp(argv[i], "-r") == 0)
		{
			++i;
			if (i > argc)
			{
				printf("Incorrect number of values for -r option\n");
				return EXIT_FAILURE;
			}

			char* str_end = nullptr;

			const float x = strtof(argv[i], &str_end);

			if (str_end == argv[i])
			{
				printf("Argument for x value of red primary could not be converted to a number\n");
				return EXIT_FAILURE;
			}

			++i;

			if (i > argc)
			{
				printf("Incorrect number of values for -r option\n");
				return EXIT_FAILURE;
			}

			const float y = strtof(argv[i], &str_end);

			if (str_end == argv[i])
			{
				printf("Argument for y value of red primary could not be converted to a number\n");
				return EXIT_FAILURE;
			}

			space.r.x = x;
			space.r.y = y;

			found[0] = true;
			
			++i;
			continue;
		}

		if (strcmp(argv[i], "-g") == 0)
		{
			++i;
			if (i > argc)
			{
				printf("Incorrect number of values for -g option\n");
				return EXIT_FAILURE;
			}

			char* str_end = nullptr;

			const float x = strtof(argv[i], &str_end);

			if (str_end == argv[i])
			{
				printf("Argument for x value of green primary could not be converted to a number\n");
				return EXIT_FAILURE;
			}

			++i;

			if (i > argc)
			{
				printf("Incorrect number of values for -g option\n");
				return EXIT_FAILURE;
			}

			const float y = strtof(argv[i], &str_end);

			if (str_end == argv[i])
			{
				printf("Argument for y value of green primary could not be converted to a number\n");
				return EXIT_FAILURE;
			}

			space.g.x = x;
			space.g.y = y;

			found[1] = true;

			++i;
			continue;
		}

		if (strcmp(argv[i], "-b") == 0)
		{
			++i;
			if (i > argc)
			{
				printf("Incorrect number of values for -b option\n");
				return EXIT_FAILURE;
			}

			char* str_end = nullptr;

			const float x = strtof(argv[i], &str_end);

			if (str_end == argv[i])
			{
				printf("Argument for x value of blue primary could not be converted to a number\n");
				return EXIT_FAILURE;
			}

			++i;

			if (i > argc)
			{
				printf("Incorrect number of values for -b option\n");
				return EXIT_FAILURE;
			}

			const float y = strtof(argv[i], &str_end);

			if (str_end == argv[i])
			{
				printf("Argument for y value of blue primary could not be converted to a number\n");
				return EXIT_FAILURE;
			}

			space.b.x = x;
			space.b.y = y;

			found[2] = true;

			++i;
			continue;
		}

		if (strcmp(argv[i], "-w") == 0)
		{
			++i;
			if (i > argc)
			{
				printf("Incorrect number of values for -w option\n");
				return EXIT_FAILURE;
			}

			char* str_end = nullptr;

			const float x = strtof(argv[i], &str_end);

			if (str_end == argv[i])
			{
				printf("Argument for x value of white point could not be converted to a number\n");
				return EXIT_FAILURE;
			}

			++i;

			if (i > argc)
			{
				printf("Incorrect number of values for -w option\n");
				return EXIT_FAILURE;
			}

			const float y = strtof(argv[i], &str_end);

			if (str_end == argv[i])
			{
				printf("Argument for y value of white point could not be converted to a number\n");
				return EXIT_FAILURE;
			}

			space.w.x = x;
			space.w.y = y;

			found[3] = true;

			++i;
			continue;
		}

		++i;
	}

	bool all_found = true;

	if (!found[0])
	{
		printf("Red primary not provided. Use option -r [x] [y] to specify red primary\n");
		all_found = false;
	}

	if (!found[1])
	{
		printf("Green primary not provided. Use option -g [x] [y] to specify green primary\n");
		all_found = false;
	}	

	if (!found[2])
	{
		printf("Blue primary not provided. Use option -b [x] [y] to specify blue primary\n");
		all_found = false;
	}

	if (!found[3])
	{
		printf("White point not provided. Use option -w [x] [y] to specify white point\n");
		all_found = false;
	}

	if (!all_found)
	{
		return EXIT_FAILURE;
	}
	
	M33 m = colour_space_to_xyz_matrix(space, verbose);

	printf("Result: \n\n");
	for (int i = 0; i < 3; ++i)
	{
		printf("%.8f, %.8f, %.8f\n", m.r[i].v[0], m.r[i].v[1], m.r[i].v[2]);
	}
}