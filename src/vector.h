#ifndef __iso3d_vector_h__
#define	__iso3d_vector_h__

/* *********************************************************************************
 * 定义三维矢量
 * define 3-dimensional vector
 * **********************************************************************************/
class Vector
{
public:
	Vector() { x = y = z = 0.; }
	Vector(double v1, double v2, double v3) {
		x = v1; y = v2; z = v3;
	}
	Vector(const Vector& v) {
		x = v.x;  y = v.y; z = v.z;
	}

	// 为满足data_io.cpp中的使用，定义新的函数
	//double X() { return this->x; }
	//double Y() { return this->y; }
	//double Z() { return this->z; }
	void setX(double x) { this->x = x; }
	void setY(double y) { this->y = y; }
	void setZ(double z) { this->z = z; }
	
	/* 矢量求模 calc. magnitude */
	double magnitude() const;
	
	/* 归一 normalization */
	void normalize();

	/* **********************************
	 * 重载操作符 overloaded operator 
	 * **********************************/

	/* 赋值 assignment */
	const Vector& operator = (const Vector& v);
	const Vector& operator = (double v);

	/* 四则运算 simple mathamatic operator */
    const Vector operator + (const Vector& v) const;
    const Vector operator - (const Vector& v) const;
    const Vector operator * (double s) const;  
    const Vector operator / (double s) const;
    
    const Vector& operator += (const Vector& v);
    const Vector& operator -= (const Vector& v);
    const Vector& operator += (double delta);
    const Vector& operator -= (double delta);

	/* 点积 & 叉积 dot & cross */
	double operator * (const Vector& v) const;
	const Vector operator ^ (const Vector& V) const;
	const Vector operator - ();

//	friend const Vector operator - (const Vector& p);
	friend const Vector operator * (double scale, const Vector& v);

public:
	double x, y, z;
};

#endif /* __iso3d_vector_h__*/
