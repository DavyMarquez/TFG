#include <mitsuba/render/shape.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/emitter.h>
#include <mitsuba/render/subsurface.h>
#include <mitsuba/render/trimesh.h>
#include <mitsuba/render/sensor.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>
#include "../dependencies/include/boost/algorithm/string.hpp"

#include "quad/QuadSampler.hpp"
#include "disk/DiskSampler.hpp"
#include "ellipse/EllipseSampler.hpp"

MTS_NAMESPACE_BEGIN //evitar problemas con otras librerias

class Ellipse : public Shape {
public:
	enum class SamplingAlgorithm {
		Area,
		Gamito,
		SolidAngle
	};

	Ellipse(const Properties &props) : Shape(props) {
		m_objectToWorld = new AnimatedTransform(props.getAnimatedTransform("toWorld", Transform()));

		if (props.getBoolean("flipNormals", false)) {
			m_objectToWorld->prependScale(Vector(1, 1, -1));
		}

		std::string sampling = boost::to_lower_copy(props.getString("sampling", "area"));
		if (sampling == "area") {
			m_sampling_type = SamplingAlgorithm::Area;
		}
		else if (sampling == "gamito") {
			m_sampling_type = SamplingAlgorithm::Gamito;
			m_sampling = EllipseSampler<double>(SphericalEllipseSampler<double>::Type::POLAR);
		}
		else if (sampling == "solidangle") {
			m_sampling_type = SamplingAlgorithm::SolidAngle;
			m_sampling = EllipseSampler<double>(SphericalEllipseSampler<double>::Type::POLAR);
		}
		else {
			Log(EError, "Unknown sampling strategy: ", sampling);
		}

		m_semiAxisA_local = props.getFloat("semiaxisa", 1.0f); //they are not scaled
		m_semiAxisB_local = props.getFloat("semiaxisb", 1.0f);
		m_invSurfaceArea = Float(0);
		m_semiAxisA = Float(0);
		m_semiAxisB = Float(0);
	}

	Ellipse(Stream *stream, InstanceManager *manager)
		:Shape(stream, manager) {
		m_objectToWorld = new AnimatedTransform(stream);
		m_sampling_type = SamplingAlgorithm::Area;

		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Shape::serialize(stream, manager);
		m_objectToWorld->serialize(stream);
	}

	void configure() {
		Shape::configure();

		const Transform &trafo = m_objectToWorld->eval(0);

		m_center = trafo(Point(0, 0, 0));
		m_semiAxisA = distance(m_center, trafo(Point(m_semiAxisA_local, 0, 0))); //in global coors
		m_semiAxisB = distance(m_center, trafo(Point(0, m_semiAxisB_local, 0)));

		m_invSurfaceArea = 1.0f / (M_PI * m_semiAxisA * m_semiAxisB);
		m_normal = normalize(trafo(Normal(0, 0, 1)));
		//cout << toString() << endl;
	}

	AABB getAABB() const {
		std::set<Float> times;
		m_objectToWorld->collectKeyframes(times);

		AABB aabb;
		for (std::set<Float>::iterator it = times.begin(); it != times.end(); ++it) {
			const Transform &trafo = m_objectToWorld->eval(*it);
			aabb.expandBy(trafo(Point(m_semiAxisA_local, 0, 0)));
			aabb.expandBy(trafo(Point(-m_semiAxisA_local, 0, 0)));
			aabb.expandBy(trafo(Point(0, m_semiAxisB_local, 0)));
			aabb.expandBy(trafo(Point(0, -m_semiAxisB_local, 0)));
		}
		return aabb;
	}

	Float getSurfaceArea() const {
		const Transform &trafo = m_objectToWorld->eval(0);
		Vector dpdu = trafo(Vector(m_semiAxisA_local, 0, 0));
		Vector dpdv = trafo(Vector(0, m_semiAxisB_local, 0));
		//CHANCES ARE I MIGHT HAVE TO CHANGE THIS 
		return M_PI * dpdu.length() * dpdv.length();
	}

	struct EllipseStorage {
		Float x, y;
	};

	inline bool rayIntersect(const Ray &_ray, Float mint, Float maxt, Float &t, void *temp) const {
		if (isHidden(_ray)) {
			return false;
		}
		Ray ray;
		m_objectToWorld->eval(ray.time).inverse().transformAffine(_ray, ray);
		Float hit = -ray.o.z / ray.d.z;

		if (!(hit >= mint && hit <= maxt))
			return false;

		Point local = ray(hit);

		if ((local.x * local.x) / (m_semiAxisA_local * m_semiAxisA_local) +
			(local.y * local.y) / (m_semiAxisB_local * m_semiAxisB_local) <= 1) {
			t = hit;

			EllipseStorage* data = static_cast<EllipseStorage*>(temp);
			if (data) {
				data->x = local.x;
				data->y = local.y;
			}
			return true;
		}
		else {
			return false;
		}
	}

	bool rayIntersect(const Ray &ray, Float mint, Float maxt) const {
		float t;
		return Ellipse::rayIntersect(ray, mint, maxt, t, NULL);
	}

	void fillIntersectionRecord(const Ray &ray,
		const void *temp, Intersection &its) const {
		const EllipseStorage* data = static_cast<const EllipseStorage*>(temp);

		Float r = std::sqrt(data->x * data->x + data->y * data->y),
			invR = (r == 0) ? 0.0f : (1.0f / r);

		Float phi = std::atan2(data->y, data->x);
		if (phi < 0)
			phi += 2 * M_PI;

		Float cosPhi = data->x * invR, sinPhi = data->y * invR;
		const Transform &trafo = m_objectToWorld->eval(ray.time);

		its.shape = this;
		if (r != 0) {
			its.dpdu = trafo(Vector(cosPhi, sinPhi, 0));
			its.dpdv = trafo(Vector(-sinPhi, cosPhi, 0));
		}
		else {
			its.dpdu = trafo(Vector(1, 0, 0));
			its.dpdv = trafo(Vector(0, 1, 0));
		}

		its.shFrame.n = normalize(trafo(Normal(0, 0, 1)));
		its.uv = Point2(r, phi * INV_TWOPI);
		its.p = ray(its.t);
		its.hasUVPartials = false;
		its.instance = nullptr;
		its.time = ray.time;
	}


	ref<TriMesh> createTriMesh() {
		const uint32_t phiSteps = 40;

		ref<TriMesh> mesh = new TriMesh(getName(),
			phiSteps - 1, 2 * phiSteps, true, true, false);

		Point *vertices = mesh->getVertexPositions();
		Normal *normals = mesh->getVertexNormals();
		Point2 *texcoords = mesh->getVertexTexcoords();
		Triangle *triangles = mesh->getTriangles();

		Float dphi = (2 * M_PI) / (Float)(phiSteps - 1);

		const Transform &trafo = m_objectToWorld->eval(0.0f);
		Point center = trafo(Point(0.0f));
		Normal normal = normalize(trafo(Normal(0, 0, 1)));

		for (uint32_t i = 0; i < phiSteps; ++i) {
			Float phi = i * dphi;
			vertices[i] = center;
			vertices[phiSteps + i] = trafo(
				Point(std::cos(phi), std::sin(phi), 0)
			);

			normals[i] = normal;
			normals[phiSteps + i] = normal;
			texcoords[i] = Point2(0.0f, phi * INV_TWOPI);
			texcoords[phiSteps + i] = Point2(1.0f, phi * INV_TWOPI);
		}

		for (uint32_t i = 0; i < phiSteps - 1; ++i) {
			triangles[i].idx[0] = i;
			triangles[i].idx[1] = i + phiSteps;
			triangles[i].idx[2] = i + phiSteps + 1;
		}

		mesh->copyAttachments(this);
		mesh->configure();

		return mesh.get();
	}

	void getNormalDerivative(const Intersection &its,
		Vector &dndu, Vector &dndv, bool shadingFrame) const {
		dndu = dndv = Vector(0.0f);
	}

	void samplePosition(PositionSamplingRecord &pRec, const Point2 &sample) const {
		const Transform &trafo = m_objectToWorld->eval(pRec.time);
		// p = [-1,1]^2
		Point2 p = warp::squareToUniformDiskConcentric(sample);
		p.x = p.x * m_semiAxisA_local; //scale disk to a ellipse
		p.y = p.y * m_semiAxisB_local;

		pRec.p = trafo(Point3(p.x, p.y, 0));

		pRec.n = normalize(trafo(Normal(0, 0, 1)));
		pRec.pdf = m_invSurfaceArea;
		pRec.measure = EArea;
	}

	Float pdfPosition(const PositionSamplingRecord &pRec) const {
		return m_invSurfaceArea;
	}

	//sample = [0, 1]^2
	void sampleDirect(DirectSamplingRecord &dRec, const Point2 &sample) const {
		//cout << "Sample Direct" << endl;
		const Transform &trafo = m_objectToWorld->eval(dRec.time);

		const Point c = trafo(Point3(0.0, 0.0, 0.0)); //center of ellipse (global coor)
		const Normal n = normalize(trafo(Normal(0.0, 0.0, 1.0))); //normal of ellipse (global coor)

		Point realRef = dRec.ref; // shading point (global coor?)

		if (m_sampling_type == SamplingAlgorithm::Area) {
			samplePosition(dRec, sample);

			// Virtual point for pdf
			dRec.d = dRec.p - dRec.ref;

			Float distSquared = dRec.d.lengthSquared();
			dRec.d /= std::sqrt(distSquared);

			Float dp = absDot(dRec.d, dRec.n);
			dRec.pdf *= (dp != 0.f) ? (distSquared / dp) : 0.f;
			dRec.measure = ESolidAngle;

			// Real point for shading
			dRec.ref = realRef;

			dRec.d = dRec.p - dRec.ref;
			dRec.dist = dRec.d.length();
			dRec.d /= dRec.dist;
		}
		else if (m_sampling_type == SamplingAlgorithm::Gamito) {
			Vector l = dRec.ref - c;
			Float d = dot(l, n);
			Point o_p = dRec.ref - d * n;

			Vector ex, ey;
			Point s;

			Vector c_d = normalize(o_p - c);
			if (c_d.lengthSquared() > Epsilon) {
				// Oriented Quad
				ey = -c_d * m_semiAxisB;
				ex = cross(n, c_d) * m_semiAxisA;

				s = c - ey - ex;
				ey = 2.f * ey;
				ex = 2.f * ex;
			}
			else {
				// Any Quad is fine
				s = trafo(Point3(-1.0, -1.0, 0.0));
				ex = trafo(Point3(1.0, -1.0, 0.0)) - s;
				ey = trafo(Point3(-1.0, 1.0, 0.0)) - s;
			}
			/*ey = trafo(Vector3(0, m_semiAxisB_local, 0));
			ex = trafo(Vector3(m_semiAxisA_local, 0, 0));
			s = trafo(Point3(-m_semiAxisA_local, -m_semiAxisB_local, 0));

			ey = 2.f * ey;
			ex = 2.f * ex;*/

			v3l::vector3<Float> o = v3l::vector3<Float>(dRec.ref.x, dRec.ref.y, dRec.ref.z);
			v3l::vector3<Float> cv = v3l::vector3<Float>(c.x, c.y, c.z);
			v3l::vector3<Float> sv = v3l::vector3<Float>(s.x, s.y, s.z);
			v3l::vector3<Float> exv = v3l::vector3<Float>(ex.x, ex.y, ex.z);
			v3l::vector3<Float> eyv = v3l::vector3<Float>(ey.x, ey.y, ey.z);

			QuadSamplingRecord<Float> samplingRecord = m_quad.createRecord(o, sv, exv, eyv);

			Float pdf = 0.0;
			v3l::vector3<Float> p = samplingRecord.sample(sample.x, sample.y, pdf);
			dRec.ref = realRef;

			v3l::vector3<Float> distance = p - cv;

			// Reject samples
			/*if ((((p.x - cv.x) * (p.x - cv.x)) / (m_semiAxisA * m_semiAxisA)) +
				(((p.y - cv.y) * (p.y - cv.y)) / (m_semiAxisB * m_semiAxisB)) > (1 + Epsilon)) {*/
			/*if(abs(p.x - cv.x) > m_semiAxisA || abs(p.y - cv.y) > m_semiAxisB){*/
			//if(length(p - cv) >= m_semiAxisA){
			//if(abs(distance.x) >= m_semiAxisA || abs(distance.y) >= m_semiAxisB){
			if (abs(distance.x) >= (m_semiAxisA + Epsilon) || abs(distance.y) >= (m_semiAxisB + Epsilon)) {
				dRec.measure = EInvalidMeasure;
				dRec.pdf = 0.0;

				return;
			}

			dRec.measure = ESolidAngle;
			dRec.pdf = this->pdfDirect(dRec);
			dRec.p = Point(p.x, p.y, p.z);
			dRec.n = n;
			dRec.d = dRec.p - dRec.ref;
			dRec.dist = dRec.d.length();
			dRec.d /= dRec.dist;

		}
		else if (m_sampling_type == SamplingAlgorithm::SolidAngle) {
			//cout << "sampling type = solid angle" << endl;
			TVector3<double> shadPoint = TVector3<double>(dRec.ref.x, dRec.ref.y, dRec.ref.z);

			//TVector3<double> shadPoint = TVector3<double>(-1, -1, -1);
			//TVector3<double> shadPoint = TVector3<double>(0, 4.9, 0);

			//cout << "shadPoint = " << shadPoint.toString() << endl;
			//Shading point transform matrix
			double toWorldMShadAux[4][4] = { {1.0, 0.0, 0.0, shadPoint.x},
											{0.0, 1.0, 0.0, shadPoint.y},
											{0.0, 0.0, 1.0, shadPoint.z},
											{0.0, 0.0, 0.0, 1.0} };

			Matrix<4, 4, double> toWorldMShad = Matrix<4, 4, double>(toWorldMShadAux);

			/*cout << "toWorldMShad" << endl;
			cout << toWorldMShad.toString() << endl;*/

			Matrix<4, 4, double> toWorldMShad_inv;
			bool useless = toWorldMShad.invert(toWorldMShad_inv);

			/*cout << "toWorldMShad_inv" << endl;
			cout << toWorldMShad_inv.toString() << endl;*/

			//Arbitrary ellipse center in global coordinates
			TVector3<double> C = TVector3<double>(m_center.x, m_center.y, m_center.z); //Ellipse center

			//cout << "C (real) = " << C.toString() << endl;

			//Arbitrary ellipse semiaxes
			double lx = m_semiAxisA,
				ly = m_semiAxisB;

			//cout << "lx = " << lx << ", ly = " << ly << endl;

			//Arbitrary ellipse trasform matrix
			double s = m_semiAxisA / m_semiAxisA_local;
			Vector n2 = trafo(Vector(0.0, 0.0, 1.0)) / s;
			TVector3<double> n = TVector3<double>(n2.x, n2.y, n2.z);
			Vector u2 = trafo(Vector(1.0, 0.0, 0.0)) / s;
			TVector3<double> u = TVector3<double>(u2.x, u2.y, u2.z);
			Vector v2 = trafo(Vector(0.0, 1.0, 0.0)) / s;
			TVector3<double> v = TVector3<double>(v2.x, v2.y, v2.z);

			double toWorldMArbAux[4][4] = { {u.x, v.x, n.x, m_center.x},
			{u.y, v.y, n.y, m_center.y},
			{u.z, v.z, n.z, m_center.z},
			{0.0, 0.0, 0.0, 1.0} };
			Matrix<4, 4, double> toWorldMArb = Matrix<4, 4, double>(toWorldMArbAux),
				toWorldMArb_real = Matrix<4, 4, double>(toWorldMArbAux);

			/*cout << "toWorldMArb_real " << endl;
			cout << toWorldMArb_real.toString() << endl;*/

			//Arbitrary ellipse transform matrix over shad point
			toWorldMArb = toWorldMShad_inv * toWorldMArb;

			/*cout << "toWorldMArb " << endl;
			cout << toWorldMArb.toString() << endl;*/

			//Arb ellipse center over shad point
			C = TVector3<double>(toWorldMArb(0, 3), toWorldMArb(1, 3), toWorldMArb(2, 3));

			//cout << "C (over shad) = " << C.toString() << endl;

			//Arb ellipse axes
			TVector3<double> Vx = TVector3<double>(toWorldMArb(0, 0), toWorldMArb(0, 1), toWorldMArb(0, 2)),
				Vy = TVector3<double>(toWorldMArb(1, 0), toWorldMArb(1, 1), toWorldMArb(1, 2)),
				Vz = TVector3<double>(toWorldMArb(2, 0), toWorldMArb(2, 1), toWorldMArb(2, 2));
			Vx = Vx / Vx.length();
			Vy = Vy / Vy.length();
			Vz = Vz / Vz.length();

			double matrixAux[4][4] = { {toWorldMArb(0,0), toWorldMArb(0,1), toWorldMArb(0,2), 0.0},
			{toWorldMArb(1,0), toWorldMArb(1,1), toWorldMArb(1,2), 0.0},
			{toWorldMArb(2,0), toWorldMArb(2,1), toWorldMArb(2,2), 0.0},
			{0.0, 0.0, 0.0, 1.0} };
			Matrix<4, 4, double> matrix = Matrix<4, 4, double>(matrixAux);

			/*cout << "matrix" << endl;
			cout << matrix.toString() << endl;*/

			double xc = C.x * Vx.x + C.y * Vx.y + C.z * Vx.z,
				yc = C.x * Vy.x + C.y * Vy.y + C.z * Vy.z,
				zc = C.x * Vz.x + C.y * Vz.y + C.z * Vz.z;

			double Qaux[3][3] = { {(zc*zc) / (lx*lx), 0, (-zc / (lx*lx))*xc},
			{0, (zc*zc) / (ly*ly), (-zc / (ly*ly))*yc},
			{ (-zc / (lx*lx))*xc, (-zc / (ly*ly))*yc,
									((xc*xc) / (lx*lx)) + ((yc*yc) / (ly*ly)) - 1} };
			Matrix<3, 3, double> Q = Matrix<3, 3, double>(Qaux);

			/*cout << "Q" << endl;
			cout << Q.toString() << endl;*/

			//eigendecomposition
			double E[3]; //eigenvalues
			Q.symEig(Q, E); //Q now stores the eigenvectors
			int order[3] = { 0, 1, 2 };
			//order Q by descending eigenvalues
			for (int i = 0; i < 3 - 1; i++) {
				int maxIndex = i;
				for (int j = i + 1; j < 3; j++) {
					if (E[j] > E[maxIndex]) {
						maxIndex = j;
					}
				}
				double aux = E[maxIndex];
				E[maxIndex] = E[i];
				E[i] = aux;
				int auxI = order[maxIndex];
				order[maxIndex] = order[i];
				order[i] = auxI;
			}
			double Qaux2[3][3] = { {Q(0,order[0]), Q(0,order[1]), Q(0,order[2])},
				{Q(1,order[0]), Q(1,order[1]), Q(1,order[2])},
				{Q(2,order[0]), Q(2,order[1]), Q(2,order[2])} };
			Q = Matrix<3, 3, double>(Qaux2);

			/*cout << "V" << endl;
			cout << Q.toString() << endl;*/

			double lambda = 1;

			//Tangent ellipse axes
			double semiAxisA_tangent = lambda * sqrt(-E[2] / E[0]),
				semiAxisB_tangent = lambda * sqrt(-E[2] / E[1]);

			/*cout << "semiAxisA_tangent = " << semiAxisA_tangent << endl;
			cout << "semiAxisB_tangent = " << semiAxisB_tangent << endl;*/

			TVector3<double> center = TVector3<double>(lambda * Q(0, 2), lambda * Q(1, 2), lambda * Q(2, 2)),
				center2 = TVector3<double>(-lambda * Q(0, 2), -lambda * Q(1, 2), -lambda * Q(2, 2));

			// "_h" refers to homogenoeus coordinates
			double centerMatrixAux[4][1] = { center[0], center[1], center[2], 1.0 };
			double centerMatrixAux2[4][1] = { center2[0], center2[1], center2[2], 1.0 };
			Matrix<4, 1, double> centerAux_h = matrix * Matrix<4, 1, double>(centerMatrixAux),
				center2Aux_h = matrix * Matrix<4, 1, double>(centerMatrixAux2);
			//cout << "centerAux_h " << centerAux_h.toString() << endl;
			TVector3<double> centerAux = TVector3<double>(centerAux_h(0, 0), centerAux_h(1, 0), centerAux_h(2, 0)),
				center2Aux = TVector3<double>(center2Aux_h(0, 0), center2Aux_h(1, 0), center2Aux_h(2, 0));

			//keep tangent ellipse center closer to arbitrary ellipse
			if ((centerAux - C).length() > (center2Aux - C).length()) {
				center = center2;
				centerAux = center2Aux;
				centerAux_h = center2Aux_h;
			}
			/*cout << "centerAux = " << centerAux_h.toString() << endl;
			cout << "center =  " << center.toString() << endl;*/

			TVector3<double> Xaxis = TVector3<double>(Q(0, 0), Q(1, 0), Q(2, 0)),
				Yaxis = TVector3<double>(Q(0, 1), Q(1, 1), Q(2, 1)),
				Zaxis = TVector3<double>(Q(0, 2), Q(1, 2), Q(2, 2));
			Xaxis = Xaxis / Xaxis.length();
			Yaxis = Yaxis / Yaxis.length();
			Zaxis = Zaxis / Yaxis.length();

			//Tangent ellipse global coordinate system 
			double Xaxis_hAux_[4][1] = { Xaxis.x, Xaxis.y, Xaxis.z, 0.0 };
			Matrix<4, 1, double> Xaxis_hAux = Matrix<4, 1, double>(Xaxis_hAux_);
			Matrix<4, 1, double> col1 = matrix * Xaxis_hAux;
			//col1 = matrix * Vector4(Xaxis.x, Xaxis.y, Xaxis.z, 0);
			//col2 = matrix * Vector4(Yaxis.x, Yaxis.y, Yaxis.z, 0);
			double Yaxis_hAux_[4][1] = { Yaxis.x, Yaxis.y, Yaxis.z, 0.0 };
			Matrix<4, 1, double> Yaxis_hAux = Matrix<4, 1, double>(Yaxis_hAux_);
			Matrix<4, 1, double> col2 = matrix * Yaxis_hAux;
			//col3 = matrix * Vector4(Zaxis.x, Zaxis.y, Zaxis.z, 0);
			double Zaxis_hAux_[4][1] = { Zaxis.x, Zaxis.y, Zaxis.z, 0.0 };
			Matrix<4, 1, double> Zaxis_hAux = Matrix<4, 1, double>(Zaxis_hAux_);
			Matrix<4, 1, double> col3 = matrix * Zaxis_hAux;
			//col4 = Vector4(centerAux.x, centerAux.y, centerAux.z, 1);
			double centerAux_hAux_[4][1] = { centerAux.x, centerAux.y, centerAux.z, 1.0 };
			Matrix<4, 1, double> centerAux_hAux = Matrix<4, 1, double>(centerAux_hAux_);
			Matrix<4, 1, double> col4 = matrix * centerAux_hAux;
			double toWorldMTangAux[4][4] = { {col1(0,0), col2(0,0), col3(0,0), centerAux.x},
			{col1(0,1), col2(0,1), col3(0,1), centerAux.y},
			{col1(0,2), col2(0,2), col3(0,2), centerAux.z},
			{col1(0,3), col2(0,3), col3(0,3), 1.0} };
			Matrix<4, 4, double> toWorldMTang = Matrix<4, 4, double>(toWorldMTangAux);
			toWorldMTang = toWorldMShad * toWorldMTang;

			/*cout << "toWorldMTang" << endl;
			cout << toWorldMTang.toString() << endl;*/

			//Tangent ellipse axes in global coordinates
			Xaxis = TVector3<double>(toWorldMTang(0, 0), toWorldMTang(1, 0), toWorldMTang(2, 0));
			Yaxis = TVector3<double>(toWorldMTang(0, 1), toWorldMTang(1, 2), toWorldMTang(2, 1));
			Zaxis = TVector3<double>(toWorldMTang(0, 2), toWorldMTang(1, 2), toWorldMTang(2, 2));

			/*Now we compute the disk that subtends the same solid angle
			  as the tangent ellipse */

			  //get pdf and if > 0, compute the disk
			dRec.measure = ESolidAngle;
			dRec.pdf = this->pdfDirect(dRec);
			//cout << "PDF = " << dRec.pdf << endl;
			//char var;
			//std::cin >> var;

			if (dRec.pdf <= 0.0) {
				dRec.pdf = 0.0;
				//cout << "return" << endl;
				dRec.measure = ESolidAngle;
				return;
			}

			double A = 1 / dRec.pdf, //spherical area
				h = A / (2 * M_PI), //height of cap
				alpha = acos(1 - A / (2 * M_PI)),
				r = sqrt(2 * h - h * h), //spherical cap radius
				a = tan(asin(r)); 


			/*cout << "A = " << A << endl;
			cout << "h = " << h << endl;
			cout << "a (radius) = " << a << endl;*/
			
			/*char var;
			std::cin >> var;*/
			//Disk radius
			double radius = a;

			//center of tangent ellipse in global coord
			centerAux = TVector3<double>(toWorldMTang(0, 3), toWorldMTang(1, 3), toWorldMTang(2, 3));

			//vector from tangent ellipse center to shad point
			TVector3<double> center2shad = shadPoint - centerAux;

			//center of disk
			TVector3<double> centerDisk = centerAux + h * center2shad;

			//cout << "centerDisk = " << centerDisk.toString() << endl;

			double toWorldMDiskAux[4][4] = { {toWorldMTang(0,0), toWorldMTang(0,1), toWorldMTang(0,2), centerDisk.x },
			{toWorldMTang(1,0), toWorldMTang(1,1), toWorldMTang(1,2), centerDisk.y},
			{toWorldMTang(2,0), toWorldMTang(2,1), toWorldMTang(2,2), centerDisk.z},
			{0.0, 0.0, 0.0, 1.0} };
			Matrix<4, 4, double> toWorldMDisk = Matrix<4, 4, double>(toWorldMDiskAux);

			/*cout << "to World M disk" << endl;
			cout << toWorldMDisk.toString() << endl;*/

			/*char var;
			std::cin >> var;*/
			//now we sample the disk
			//const Vector refToCenter = centerDisk - dRec.ref;
			const Vector refToCenter = Point(centerDisk.x, centerDisk.y, centerDisk.z) -
				Point(shadPoint.x, shadPoint.y, shadPoint.z);
			const double refDist2 = refToCenter.lengthSquared();
			const double invRefDist = static_cast<double>(1) / std::sqrt(refDist2);

			double cosAlpha = cos(alpha);

			dRec.d = Frame(refToCenter * invRefDist).toWorld(
				warp::squareToUniformCone(cosAlpha, sample));

			dRec.d = dRec.d / dRec.d.length();

			//cout << "sampled vector" << endl;
			//cout << dRec.d.x << ", " << dRec.d.y << ", " << dRec.d.z << endl;

			//check collision with disk plane
			TVector3<double> planeNormal = TVector3<double>(toWorldMDisk(0, 2), toWorldMDisk(1, 2), toWorldMDisk(2, 2));
			planeNormal = planeNormal / planeNormal.length();

			//cout << "planeNormal" << endl;
			/*cout << planeNormal.toString() << endl;
			cout << "centerDisk" << endl;
			cout << centerDisk.toString() << endl;*/
			TVector3<double> l = TVector3<double>(dRec.d.x, dRec.d.y, dRec.d.z);
			l = l / l.length();
			//cout << l.toString() << endl;
			if (dot(l, planeNormal) < 0) {
				planeNormal = -planeNormal;
				//cout << "change normal direction" << endl;
			}
			//cout << planeNormal.toString() << endl;
			double denom = dot(planeNormal, l),
				t = 0;
			if (denom > 0) {
				TVector3<double> shad2disk = centerDisk - shadPoint;
				/*cout << "shad2disk" << endl;
				cout << shad2disk.toString() << endl;*/
				t = dot(shad2disk, planeNormal) / denom; //dist to disk plane
				if (t < 0) {
					dRec.measure = ESolidAngle;
					dRec.pdf = 0.0;
					//cout << "t < 0" << endl;
					/*char variable;
					std::cin >> variable;*/
					return;
				}
			}
			else {
				//cout << "denom = 0" << endl;
				dRec.measure = ESolidAngle;
				dRec.pdf = 0.0;
				/*char variable;
				std::cin >> variable;*/

				return;
			}

			/*if we are here, that means it collided the disk plane, now check
			  if it hits the disk */
			  //cout << "collided with plane" << endl;
			  //Point hit = dRec.ref + l * t;
			TVector3<double> hit = shadPoint + l * t;
			TVector3<double> hit2center = hit - centerDisk;
			if (hit2center.length() > radius) {
				dRec.measure = ESolidAngle;
				dRec.pdf = 0.0;
				//cout << "not collision with disk" << endl;
				/*char variable;
				std::cin >> variable;*/
				return;
			}
			//cout << "collided with disk" << endl;
			TVector3<double> hitAux = hit;
			/*hit = Point(centerDisk.x + radius , centerDisk.y, centerDisk.z);*/
			//cout << "disk hit: " << hit.toString() << endl;

			//now map disk point to tangent ellipse point
			Matrix<4, 4, double> toWorldMDisk_inv;
			useless = toWorldMDisk.invert(toWorldMDisk_inv);

			double hit_h_aux[4] = { hit.x, hit.y, hit.z, 1 };
			Matrix<4, 1, double> hit_h = Matrix<4, 1, double>(hit_h_aux);
			hit_h = toWorldMDisk_inv * hit_h;
			double hit_h_aux2[4] = { (hit_h(0,0) / radius) * semiAxisA_tangent,
				(hit_h(1,0) / radius) * semiAxisB_tangent, hit_h(2,0), 1.0 };
			hit_h = Matrix<4, 1, double>(hit_h_aux2);

			//hit on the tangent ellipse
			hit_h = toWorldMTang * hit_h;
			hit = TVector3<double>(hit_h(0, 0), hit_h(1, 0), hit_h(2, 0));

			//cout << "tangent hit = " << hit.toString() << endl;

			//now collide with arbitrary ellipse
			l = (hit - shadPoint) / (hit - shadPoint).length();

			TVector3<double> ellipseNormal = TVector3<double>(toWorldMArb_real(0, 2),
				toWorldMArb_real(1, 2),
				toWorldMArb_real(2, 2));
			ellipseNormal = ellipseNormal / ellipseNormal.length();
			ellipseNormal = -ellipseNormal;

			/*cout << "ellipse normal" << endl;
			cout << ellipseNormal.toString() << endl;*/
			/*if (dot(l, ellipseNormal) < 0.0) {
				ellipseNormal = -ellipseNormal;
			}*/
			//ellipseNormal = -ellipseNormal;
			/*cout << "ellipse normal negated" << endl;
			cout << ellipseNormal.toString() << endl;*/
			denom = dot(ellipseNormal, l);
			//cout << "denom = " << denom << endl;
			TVector3<double> ellipseCenter = TVector3<double>(m_center.x, m_center.y, m_center.z);
			if (denom > 0) {
				TVector3<double> p0l0 = (ellipseCenter - hit);
				/*cout << "p0l0" << endl;
				cout << p0l0.toString() << endl;*/
				double t = dot(p0l0, ellipseNormal) / denom;
				//cout << "t = " << t << endl;
				if (t != 0) {
					//cout << "choca con la elipse" << endl;
					hit = hit + t * l;
				}
				else {
					//cout << "no choca con la elipse arbitraria" << endl;
					dRec.measure = ESolidAngle;
					dRec.pdf = 0.0;
					return;
				}
			}
			else {
				//cout << "denom < 0" << endl;
				dRec.measure = ESolidAngle;
				dRec.pdf = 0.0;		
				return;
			}
			//cout << "arb hit = " << hit.toString() << endl;
			dRec.measure = ESolidAngle;
			dRec.pdf = this->pdfDirect(dRec);
			if (dRec.pdf <= 0.0) {
				dRec.measure = ESolidAngle;
				dRec.pdf = 0.0;
				return;
			}
			char var;
			//std::cin >> var;
			//cout << "shading point = " << shadPoint.toString() << endl;
			dRec.dist = (hit - shadPoint).length();
			dRec.n = normalize(trafo(Normal(0.0, 0.0, 1.0)));
			dRec.p = Point(hit.x, hit.y, hit.z);
			dRec.d = dRec.p - dRec.ref;
			dRec.d /= dRec.dist; 
			dRec.ref = realRef;
			//cout << dRec.p.toString() << endl;




		}

		return;
	}

	Float pdfDirect(const DirectSamplingRecord &dRec) const {
		//cout << "pdf Direct" << endl;
		const Transform &trafo = m_objectToWorld->eval(dRec.time);

		const Point c = trafo(Point3(0.0, 0.0, 0.0));
		const Normal n = trafo(Normal(0.0, 0.0, 1.0));
		const Normal u = trafo(Normal(m_semiAxisA_local, 0.0, 0.0));
		const Normal v = trafo(Normal(0.0, m_semiAxisB_local, 0.0));

		DirectSamplingRecord lRec(dRec);
		Float d = dot((dRec.ref - c), n);

		// Use the default procedure
		if (m_sampling_type == SamplingAlgorithm::Area) {
			Float pdfPos = pdfPosition(lRec);

			if (lRec.measure == ESolidAngle) {
				Vector d = lRec.p - lRec.ref;

				Float distSquared = d.lengthSquared();
				d /= std::sqrt(distSquared);

				Float dp = absDot(d, lRec.n);
				return pdfPos * ((dp != 0.f) ? (distSquared / dp) : 0.f);
			}
			else if (lRec.measure == EArea) {
				return pdfPos;
			}
		}
		else if (m_sampling_type == SamplingAlgorithm::Gamito ||
			m_sampling_type == SamplingAlgorithm::SolidAngle) {
			if (lRec.measure == ESolidAngle) {
				// In order to apply the Heitz transform we need the 
				// coordinate system without normalization
				Float s = m_semiAxisA / m_semiAxisA_local;
				Vector n = trafo(Vector(0.0, 0.0, 1.0)) / s;
				Vector u = trafo(Vector(1.0, 0.0, 0.0)) / s;
				Vector v = trafo(Vector(0.0, 1.0, 0.0)) / s;
				//shading point
				v3l::vector3<double> o = v3l::vector3<double>(lRec.ref.x, lRec.ref.y, lRec.ref.z);
				//v3l::vector3<double> o = v3l::vector3<double>(0, 4.9, 0);

				//ellipse center
				v3l::vector3<double> cd = v3l::vector3<double>(c.x, c.y, c.z);
				//ellipse normal
				v3l::vector3<double> nd = v3l::vector3<double>(n.x, n.y, n.z);
				//ellipse x axis 
				v3l::vector3<double> ud = v3l::vector3<double>(u.x, u.y, u.z);
				//ellipse y axis 
				v3l::vector3<double> vd = v3l::vector3<double>(v.x, v.y, v.z);


				EllipseSamplingRecord<double> samplingRecord = m_sampling.createRecord(o,
					m_semiAxisA, m_semiAxisB,
					ud, vd, nd, cd);

				v3l::vector3<double> p = v3l::vector3<double>(lRec.p.x, lRec.p.y, lRec.p.z);
				Float pdf = samplingRecord.samplePdf(p);
				//cout << p.x << ", " << p.y << ", " <<p.z << endl;
				//cout << "pdf: " << pdf << endl;


				return pdf;
			}
			else if (lRec.measure == EArea) {
				return pdfPosition(lRec);
			}
		}

		return 0.0;
	}

	size_t getPrimitiveCount() const {
		return 1;
	}

	size_t getEffectivePrimitiveCount() const {
		return 1;

	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "Ellipse[" << endl
			<< " objectToWorld = " << indent(m_objectToWorld->toString()) << ", " << endl
			<< " center = {" << m_center.x << ", " << m_center.y << ", " << m_center.z << "}, " << endl
			<< " semiAxisA = " << m_semiAxisA << ", " << endl
			<< " semiAxisB = " << m_semiAxisB << endl
			<< " invSurfaceArea = " << m_invSurfaceArea << endl
			<< " bsdf = " << indent(m_bsdf.toString()) << ", " << endl;
		if (isMediumTransition()) {
			oss << "  interiorMedium = " << indent(m_interiorMedium.toString()) << "," << endl
				<< "  exteriorMedium = " << indent(m_exteriorMedium.toString()) << "," << endl;
		}
		oss << " emitter = " << indent(m_emitter.toString()) << "," << endl
			<< " sensor = " << indent(m_sensor.toString()) << "," << endl
			<< " subsurface = " << indent(m_subsurface.toString()) << endl
			<< "]" << endl;
		return oss.str();
	}

	MTS_DECLARE_CLASS() //para que mitsuba la reconozca como clase nativa
private:
	Point m_center; //Ellipse center (global)
	Normal m_normal; //Ellipse normal
	Float m_semiAxisA; // Mayor Semiaxis A (global)
	Float m_semiAxisB; //Minor Semiaxis B (global)
	Float m_semiAxisA_local; // Mayor Semiaxis A (local)
	Float m_semiAxisB_local; //Minor Semiaxis B (local)
	Float m_invSurfaceArea;

	ref<AnimatedTransform> m_objectToWorld; //transformation matrix
	SamplingAlgorithm m_sampling_type;
	EllipseSampler<double> m_sampling;
	QuadSampler<Float> m_quad;
};

MTS_IMPLEMENT_CLASS(Ellipse, false, Shape) //para que mitsuba la reconozca como clase nativa
MTS_EXPORT_PLUGIN(Ellipse, "Ellipse intersection primitive")
MTS_NAMESPACE_END
