//
// File: TransformConstraintNode.cpp
//
// Dependency Graph Node: transformConstraint
//
// Author: Ben Singleton
//

#include "TransformConstraint.h"


MObject TransformConstraint::restTranslate;
MObject TransformConstraint::restTranslateX;
MObject TransformConstraint::restTranslateY;
MObject TransformConstraint::restTranslateZ;
MObject TransformConstraint::restRotate;
MObject TransformConstraint::restRotateX;
MObject TransformConstraint::restRotateY;
MObject TransformConstraint::restRotateZ;
MObject TransformConstraint::restScale;
MObject TransformConstraint::restScaleX;
MObject TransformConstraint::restScaleY;
MObject TransformConstraint::restScaleZ;

MObject TransformConstraint::target;
MObject TransformConstraint::targetWeight;
MObject TransformConstraint::targetParentMatrix;
MObject TransformConstraint::targetTranslate;
MObject TransformConstraint::targetTranslateX;
MObject TransformConstraint::targetTranslateY;
MObject TransformConstraint::targetTranslateZ;
MObject TransformConstraint::targetOffsetTranslate;
MObject TransformConstraint::targetOffsetTranslateX;
MObject TransformConstraint::targetOffsetTranslateY;
MObject TransformConstraint::targetOffsetTranslateZ;
MObject TransformConstraint::targetJointOrient;
MObject TransformConstraint::targetJointOrientX;
MObject TransformConstraint::targetJointOrientY;
MObject TransformConstraint::targetJointOrientZ;
MObject TransformConstraint::targetRotate;
MObject TransformConstraint::targetRotateX;
MObject TransformConstraint::targetRotateY;
MObject TransformConstraint::targetRotateZ;
MObject TransformConstraint::targetOffsetRotate;
MObject TransformConstraint::targetOffsetRotateX;
MObject TransformConstraint::targetOffsetRotateY;
MObject TransformConstraint::targetOffsetRotateZ;
MObject TransformConstraint::targetRotateOrder;
MObject TransformConstraint::targetScale;
MObject TransformConstraint::targetScaleX;
MObject TransformConstraint::targetScaleY;
MObject TransformConstraint::targetScaleZ;
MObject TransformConstraint::targetOffsetScale;
MObject TransformConstraint::targetOffsetScaleX;
MObject TransformConstraint::targetOffsetScaleY;
MObject TransformConstraint::targetOffsetScaleZ;
MObject TransformConstraint::targetRotatePivot;
MObject TransformConstraint::targetRotatePivotX;
MObject TransformConstraint::targetRotatePivotY;
MObject TransformConstraint::targetRotatePivotZ;
MObject TransformConstraint::targetRotatePivotTranslate;
MObject TransformConstraint::targetRotatePivotTranslateX;
MObject TransformConstraint::targetRotatePivotTranslateY;
MObject TransformConstraint::targetRotatePivotTranslateZ;
MObject TransformConstraint::targetScalePivot;
MObject TransformConstraint::targetScalePivotX;
MObject TransformConstraint::targetScalePivotY;
MObject TransformConstraint::targetScalePivotZ;
MObject TransformConstraint::targetScalePivotTranslate;
MObject TransformConstraint::targetScalePivotTranslateX;
MObject TransformConstraint::targetScalePivotTranslateY;
MObject TransformConstraint::targetScalePivotTranslateZ;

MObject TransformConstraint::constraintTranslate;
MObject TransformConstraint::constraintTranslateX;
MObject TransformConstraint::constraintTranslateY;
MObject TransformConstraint::constraintTranslateZ;
MObject TransformConstraint::constraintRotate;
MObject TransformConstraint::constraintRotateX;
MObject TransformConstraint::constraintRotateY;
MObject TransformConstraint::constraintRotateZ;
MObject TransformConstraint::constraintRotateOrder;
MObject TransformConstraint::constraintJointOrient;
MObject TransformConstraint::constraintJointOrientX;
MObject TransformConstraint::constraintJointOrientY;
MObject TransformConstraint::constraintJointOrientZ;
MObject TransformConstraint::constraintScale;
MObject TransformConstraint::constraintScaleX;
MObject TransformConstraint::constraintScaleY;
MObject TransformConstraint::constraintScaleZ;
MObject TransformConstraint::constraintMatrix;
MObject TransformConstraint::constraintInverseMatrix;
MObject TransformConstraint::constraintWorldMatrix;
MObject TransformConstraint::constraintWorldInverseMatrix;
MObject TransformConstraint::constraintParentInverseMatrix;

MString	TransformConstraint::inputCategory("Input");
MString	TransformConstraint::restCategory("Rest");
MString	TransformConstraint::targetCategory("Target");
MString	TransformConstraint::outputCategory("Output");

MString TransformConstraint::classification("animation");
MTypeId	TransformConstraint::id(0x0013b1c2);


template<class N> N lerp(N start, N end, double weight)
/**
Linearly interpolates the two given numbers using the supplied weight.

@param start: The start number.
@param end: The end number.
@param weight: The amount to blend.
@return: The interpolated value.
*/
{

	return (start * (1.0 - weight)) + (end * weight);

};


TransformConstraint::TransformConstraint() {};
TransformConstraint::~TransformConstraint() {};


MStatus TransformConstraint::compute(const MPlug& plug, MDataBlock& data)
/**
This method should be overridden in user defined nodes.
Recompute the given output based on the nodes inputs.
The plug represents the data value that needs to be recomputed, and the data block holds the storage for all of the node's attributes.
The MDataBlock will provide smart handles for reading and writing this node's attribute values.
Only these values should be used when performing computations!

@param plug: Plug representing the attribute that needs to be recomputed.
@param data: Data block containing storage for the node's attributes.
@return: Return status.
*/
{

	MStatus status;

	// Check requested attribute
	//
	MObject attribute = plug.attribute(&status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	MFnAttribute fnAttribute(attribute, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	if (fnAttribute.hasCategory(TransformConstraint::outputCategory))
	{

		// Get input data handles
		//
		MDataHandle restTranslateHandle = data.inputValue(TransformConstraint::restTranslate, &status);
		CHECK_MSTATUS_AND_RETURN_IT(status);

		MDataHandle restRotateHandle = data.inputValue(TransformConstraint::restRotate, &status);
		CHECK_MSTATUS_AND_RETURN_IT(status);

		MDataHandle restScaleHandle = data.inputValue(TransformConstraint::restScale, &status);
		CHECK_MSTATUS_AND_RETURN_IT(status);

		MDataHandle constraintParentInverseMatrixHandle = data.inputValue(TransformConstraint::constraintParentInverseMatrix, &status);
		CHECK_MSTATUS_AND_RETURN_IT(status);

		MDataHandle constraintRotateOrderHandle = data.inputValue(TransformConstraint::constraintRotateOrder, &status);
		CHECK_MSTATUS_AND_RETURN_IT(status);

		MDataHandle constraintJointOrientHandle = data.inputValue(TransformConstraint::constraintJointOrient, &status);
		CHECK_MSTATUS_AND_RETURN_IT(status);

		MArrayDataHandle targetArrayHandle = data.inputArrayValue(TransformConstraint::target, &status);
		CHECK_MSTATUS_AND_RETURN_IT(status);

		// Get values from handles
		//
		MEulerRotation::RotationOrder constraintRotateOrder = MEulerRotation::RotationOrder(constraintRotateOrderHandle.asShort());

		MEulerRotation constraintJointOrient = MEulerRotation(constraintJointOrientHandle.asVector(), constraintRotateOrder);
		MMatrix constraintJointOrientMatrix = constraintJointOrient.asMatrix();

		MMatrix constraintParentInverseMatrix = constraintParentInverseMatrixHandle.asMatrix();
		MMatrix constraintParentMatrix = constraintParentInverseMatrix.inverse();

		MVector restTranslate = restTranslateHandle.asVector();
		MMatrix restTranslateMatrix = TransformConstraint::createPositionMatrix(restTranslate);
		MEulerRotation restRotate = MEulerRotation(restRotateHandle.asVector(), constraintRotateOrder);
		MMatrix restRotateMatrix = restRotate.asMatrix();
		MVector restScale = restScaleHandle.asVector();
		MMatrix restScaleMatrix = TransformConstraint::createScaleMatrix(restScale);
		MMatrix restMatrix = restScaleMatrix * restRotateMatrix * restTranslateMatrix;
		MMatrix restWorldMatrix = restMatrix * constraintParentInverseMatrix.inverse();

		// Collect target matrices
		//
		unsigned int targetCount = targetArrayHandle.elementCount();

		MFloatArray targetWeights = MFloatArray(targetCount);
		MMatrixArray targetMatrices = MMatrixArray(targetCount);
		MMatrixArray targetWorldMatrices = MMatrixArray(targetCount);

		MDataHandle targetHandle;
		MDataHandle targetWeightHandle, targetParentMatrixHandle;
		MDataHandle targetTranslateHandle, targetOffsetTranslateHandle;
		MDataHandle targetJointOrientHandle, targetRotateHandle, targetOffsetRotateHandle, targetRotateOrderHandle;
		MDataHandle targetScaleHandle, targetOffsetScaleHandle;
		MDataHandle targetRotatePivotHandle, targetRotatePivotTranslateHandle;
		MDataHandle targetScalePivotHandle, targetScalePivotTranslateHandle;

		MEulerRotation::RotationOrder targetRotateOrder;
		MMatrix targetOffsetMatrix;
		MMatrix targetParentMatrix;
		MMatrix targetTranslateMatrix, targetOffsetTranslateMatrix;
		MMatrix targetJointOrientMatrix, targetRotateMatrix, targetOffsetRotateMatrix;
		MMatrix targetScaleMatrix, targetOffsetScaleMatrix;
		MMatrix targetRotatePivotMatrix, targetRotatePivotTranslateMatrix;
		MMatrix	targetScalePivotMatrix, targetScalePivotTranslateMatrix;

		for (unsigned int i = 0; i < targetCount; i++)
		{

			// Jump to array element
			//
			status = targetArrayHandle.jumpToArrayElement(i);
			CHECK_MSTATUS_AND_RETURN_IT(status)

			targetHandle = targetArrayHandle.inputValue(&status);
			CHECK_MSTATUS_AND_RETURN_IT(status);

			// Get target data handles
			//
			targetWeightHandle = targetHandle.child(TransformConstraint::targetWeight);
			targetParentMatrixHandle = targetHandle.child(TransformConstraint::targetParentMatrix);
			targetTranslateHandle = targetHandle.child(TransformConstraint::targetTranslate);
			targetOffsetTranslateHandle = targetHandle.child(TransformConstraint::targetOffsetTranslate);
			targetJointOrientHandle = targetHandle.child(TransformConstraint::targetJointOrient);
			targetRotateHandle = targetHandle.child(TransformConstraint::targetRotate);
			targetOffsetRotateHandle = targetHandle.child(TransformConstraint::targetOffsetRotate);
			targetRotateOrderHandle = targetHandle.child(TransformConstraint::targetRotateOrder);
			targetScaleHandle = targetHandle.child(TransformConstraint::targetScale);
			targetOffsetScaleHandle = targetHandle.child(TransformConstraint::targetOffsetScale);
			targetRotatePivotHandle = targetHandle.child(TransformConstraint::targetRotatePivot);
			targetRotatePivotTranslateHandle = targetHandle.child(TransformConstraint::targetRotatePivotTranslate);
			targetScalePivotHandle = targetHandle.child(TransformConstraint::targetScalePivot);
			targetScalePivotTranslateHandle = targetHandle.child(TransformConstraint::targetScalePivotTranslate);

			// Get weight value
			//
			targetWeights[i] = targetWeightHandle.asFloat();

			// Get rotate order
			//
			targetRotateOrder = MEulerRotation::RotationOrder(targetRotateOrderHandle.asShort());

			// Compute target offset matrix
			//
			targetOffsetTranslateMatrix = TransformConstraint::createPositionMatrix(targetOffsetTranslateHandle.asVector());
			targetOffsetRotateMatrix = MEulerRotation(targetOffsetRotateHandle.asVector(), targetRotateOrder).asMatrix();
			targetOffsetScaleMatrix = TransformConstraint::createScaleMatrix(targetOffsetScaleHandle.asVector());

			targetOffsetMatrix = targetOffsetScaleMatrix * targetOffsetRotateMatrix * targetOffsetTranslateMatrix;

			// Compute target matrices
			//
			targetParentMatrix = targetParentMatrixHandle.asMatrix();
			targetTranslateMatrix = TransformConstraint::createPositionMatrix(targetTranslateHandle.asVector());
			
			targetJointOrientMatrix = MEulerRotation(targetJointOrientHandle.asVector(), targetRotateOrder).asMatrix();
			targetRotateMatrix = MEulerRotation(targetRotateHandle.asVector(), targetRotateOrder).asMatrix();
			
			targetScaleMatrix = TransformConstraint::createScaleMatrix(targetScaleHandle.asDouble3());
			
			targetRotatePivotMatrix = TransformConstraint::createPositionMatrix(targetRotatePivotHandle.asVector());
			targetRotatePivotTranslateMatrix = TransformConstraint::createPositionMatrix(targetRotatePivotTranslateHandle.asVector());
			targetScalePivotMatrix = TransformConstraint::createPositionMatrix(targetScalePivotHandle.asVector());
			targetScalePivotTranslateMatrix = TransformConstraint::createPositionMatrix(targetScalePivotTranslateHandle.asVector());

			targetMatrices[i] = targetScalePivotMatrix.inverse() * targetScaleMatrix * targetScalePivotMatrix * targetScalePivotTranslateMatrix * targetRotatePivotMatrix.inverse() * targetJointOrientMatrix * targetRotateMatrix * targetRotatePivotMatrix * targetRotatePivotTranslateMatrix * targetTranslateMatrix;
			targetWorldMatrices[i] = (targetOffsetMatrix * targetMatrices[i]) * targetParentMatrix;

		}

		// Get output data handles
		//
		MDataHandle constraintTranslateXHandle = data.outputValue(TransformConstraint::constraintTranslateX, &status);
		CHECK_MSTATUS_AND_RETURN_IT(status);

		MDataHandle constraintTranslateYHandle = data.outputValue(TransformConstraint::constraintTranslateY, &status);
		CHECK_MSTATUS_AND_RETURN_IT(status);

		MDataHandle constraintTranslateZHandle = data.outputValue(TransformConstraint::constraintTranslateZ, &status);
		CHECK_MSTATUS_AND_RETURN_IT(status);

		MDataHandle constraintRotateXHandle = data.outputValue(TransformConstraint::constraintRotateX, &status);
		CHECK_MSTATUS_AND_RETURN_IT(status);

		MDataHandle constraintRotateYHandle = data.outputValue(TransformConstraint::constraintRotateY, &status);
		CHECK_MSTATUS_AND_RETURN_IT(status);

		MDataHandle constraintRotateZHandle = data.outputValue(TransformConstraint::constraintRotateZ, &status);
		CHECK_MSTATUS_AND_RETURN_IT(status);

		MDataHandle constraintScaleXHandle = data.outputValue(TransformConstraint::constraintScaleX, &status);
		CHECK_MSTATUS_AND_RETURN_IT(status);

		MDataHandle constraintScaleYHandle = data.outputValue(TransformConstraint::constraintScaleY, &status);
		CHECK_MSTATUS_AND_RETURN_IT(status);

		MDataHandle constraintScaleZHandle = data.outputValue(TransformConstraint::constraintScaleZ, &status);
		CHECK_MSTATUS_AND_RETURN_IT(status);

		MDataHandle constraintMatrixHandle = data.outputValue(TransformConstraint::constraintMatrix, &status);
		CHECK_MSTATUS_AND_RETURN_IT(status);

		MDataHandle constraintInverseMatrixHandle = data.outputValue(TransformConstraint::constraintInverseMatrix, &status);
		CHECK_MSTATUS_AND_RETURN_IT(status);

		MDataHandle constraintWorldMatrixHandle = data.outputValue(TransformConstraint::constraintWorldMatrix, &status);
		CHECK_MSTATUS_AND_RETURN_IT(status);

		MDataHandle constraintWorldInverseMatrixHandle = data.outputValue(TransformConstraint::constraintWorldInverseMatrix, &status);
		CHECK_MSTATUS_AND_RETURN_IT(status);

		// Compute weighted constraint matrix
		//
		MMatrix constraintWorldMatrix = TransformConstraint::blendMatrices(restWorldMatrix, targetWorldMatrices, targetWeights);
		MMatrix constraintMatrix = constraintWorldMatrix * constraintParentInverseMatrix;

		// Set translation constraint
		//
		MVector constraintTranslate = TransformConstraint::getTranslationPart(constraintMatrix);

		constraintTranslateXHandle.setMDistance(MDistance(constraintTranslate.x, MDistance::kCentimeters));
		constraintTranslateYHandle.setMDistance(MDistance(constraintTranslate.y, MDistance::kCentimeters));
		constraintTranslateZHandle.setMDistance(MDistance(constraintTranslate.z, MDistance::kCentimeters));

		constraintTranslateXHandle.setClean();
		constraintTranslateYHandle.setClean();
		constraintTranslateZHandle.setClean();

		// Set rotation constraint
		//
		MMatrix constraintRotateMatrix = TransformConstraint::createRotationMatrix(constraintMatrix) * constraintJointOrientMatrix.inverse();

		MEulerRotation constraintRotate = TransformConstraint::getRotationPart(constraintRotateMatrix).asEulerRotation();
		constraintRotate.reorderIt(constraintRotateOrder);

		constraintRotateXHandle.setMAngle(MAngle(constraintRotate.x, MAngle::kRadians));
		constraintRotateYHandle.setMAngle(MAngle(constraintRotate.y, MAngle::kRadians));
		constraintRotateZHandle.setMAngle(MAngle(constraintRotate.z, MAngle::kRadians));

		constraintRotateXHandle.setClean();
		constraintRotateYHandle.setClean();
		constraintRotateZHandle.setClean();

		// Set scale constraint
		//
		MVector constraintScale = TransformConstraint::getScalePart(constraintMatrix);

		constraintScaleXHandle.setDouble(constraintScale.x);
		constraintScaleYHandle.setDouble(constraintScale.y);
		constraintScaleZHandle.setDouble(constraintScale.z);

		constraintScaleXHandle.setClean();
		constraintScaleYHandle.setClean();
		constraintScaleZHandle.setClean();

		// Commit matrices to handles
		//
		constraintMatrixHandle.setMMatrix(constraintMatrix);
		constraintInverseMatrixHandle.setMMatrix(constraintMatrix.inverse());
		constraintWorldMatrixHandle.setMMatrix(constraintWorldMatrix);
		constraintWorldInverseMatrixHandle.setMMatrix(constraintWorldMatrix.inverse());

		constraintMatrixHandle.setClean();
		constraintInverseMatrixHandle.setClean();
		constraintWorldMatrixHandle.setClean();
		constraintWorldInverseMatrixHandle.setClean();

		// Mark data block as clean
		//
		status = data.setClean(plug);
		CHECK_MSTATUS_AND_RETURN_IT(status);

		return MS::kSuccess;

	}
	else
	{
		
		return MS::kUnknownParameter;

	}

}


MStatus TransformConstraint::decomposeTransformMatrix(const MMatrix& matrix, MVector& position, MQuaternion& rotation, MVector& scale)
/**
Returns the translate, rotate and scale components from the supplied transform matrix.

@param matrix: The transform matrix to extract from.
@param position: The translation component.
@param rotation: The rotation component.
@param scale: The scale component.
@return: Status code.
*/
{

	position = TransformConstraint::getTranslationPart(matrix);
	rotation = TransformConstraint::getRotationPart(matrix);
	scale = TransformConstraint::getScalePart(matrix);

	return MS::kSuccess;

};


MVector TransformConstraint::getTranslationPart(const MMatrix& matrix)
/**
Returns the translation component from the supplied transform matrix.

@param matrix: The transform matrix to extract from.
@return: The translation component.
*/
{

	return MVector(matrix(3, 0), matrix(3, 1), matrix(3, 2));

};


MQuaternion TransformConstraint::getRotationPart(const MMatrix& matrix)
/**
Returns the rotation component from the supplied transform matrix.

@param matrix: The transform matrix to extract from.
@return: The rotation component.
*/
{

	MQuaternion rotationPart;
	rotationPart = TransformConstraint::createRotationMatrix(matrix);

	return rotationPart;

};


MVector TransformConstraint::getScalePart(const MMatrix& matrix)
/**
Returns the scale component from the supplied transform matrix.

@param matrix: The transform matrix to extract from.
@return: The scale component.
*/
{

	MVector xAxis = MVector(matrix(0, 0), matrix(0, 1), matrix(0, 2));
	MVector yAxis = MVector(matrix(1, 0), matrix(1, 1), matrix(1, 2));
	MVector zAxis = MVector(matrix(2, 0), matrix(2, 1), matrix(2, 2));

	return MVector(xAxis.length(), yAxis.length(), zAxis.length());

};


MMatrix TransformConstraint::createPositionMatrix(const MPoint& position)
/**
Creates a position matrix from the given point.

@param position: The point to convert.
@return: The new position matrix.
*/
{

	double matrixRows[4][4] = {
		{ 1.0, 0.0, 0.0, 0.0 },
		{ 0.0, 1.0, 0.0, 0.0 },
		{ 0.0, 0.0, 1.0, 0.0 },
		{ position.x, position.y, position.z, position.w },
	};

	return MMatrix(matrixRows);

};


MMatrix TransformConstraint::createPositionMatrix(const MVector& position)
/**
Creates a position matrix from the given vector.

@param position: The vector to convert.
@return: The new position matrix.
*/
{

	return TransformConstraint::createPositionMatrix(MPoint(position));

};


MMatrix TransformConstraint::createPositionMatrix(const MMatrix& matrix)
/**
Returns the position component from the supplied transform matrix.

@param position: The transform matrix to extract from.
@return: The new position matrix.
*/
{

	return TransformConstraint::createPositionMatrix(TransformConstraint::getTranslationPart(matrix));

};


MMatrix TransformConstraint::createRotationMatrix(const MMatrix& matrix)
/**
Returns the rotation component from the supplied transform matrix.

@param matrix: The transform matrix to extract from.
@return: The new rotation matrix.
*/
{

	MVector xAxis = MVector(matrix(0, 0), matrix(0, 1), matrix(0, 2)).normal();
	MVector yAxis = MVector(matrix(1, 0), matrix(1, 1), matrix(1, 2)).normal();
	MVector zAxis = MVector(matrix(2, 0), matrix(2, 1), matrix(2, 2)).normal();

	double matrixRows[4][4] = {
		{ xAxis.x, xAxis.y, xAxis.z, 0.0 },
		{ yAxis.x, yAxis.y, yAxis.z, 0.0 },
		{ zAxis.x, zAxis.y, zAxis.z, 0.0 },
		{ 0.0, 0.0, 0.0, 1.0 },
	};

	return MMatrix(matrixRows);

};


MMatrix TransformConstraint::createScaleMatrix(const MVector& scale)
/**
Returns a scale matrix from the supplied vector.

@param scale: The vector to convert.
@return: The new scale matrix.
*/
{

	double matrixRows[4][4] = {
		{ scale.x, 0.0, 0.0, 0.0 },
		{ 0.0, scale.y, 0.0, 0.0 },
		{ 0.0, 0.0, scale.z, 0.0 },
		{ 0.0, 0.0, 0.0, 1.0 },
	};

	return MMatrix(matrixRows);

};


MMatrix TransformConstraint::createScaleMatrix(const MMatrix& matrix)
/**
Returns a scale matrix from the supplied transform matrix.

@param matrix: The transform matrix to extract from.
@return: The new scale matrix.
*/
{

	return TransformConstraint::createScaleMatrix(TransformConstraint::getScalePart(matrix));

};


MMatrix TransformConstraint::blendMatrices(const MMatrix& startMatrix, const MMatrix& endMatrix, const float weight)
/**
Interpolates the two given matrices using the supplied weight.
Both translate and scale will be lerp'd while rotation will be slerp'd.

@param startMatrix: The start matrix.
@param endMatrix: The end matrix.
@param weight: The amount to blend.
@return: The interpolated matrix.
*/
{

	MStatus status;

	// Decompose transform matrices
	//
	MVector startTranslation, endTranslation;
	MVector startScale, endScale;
	MQuaternion startQuat, endQuat;

	TransformConstraint::decomposeTransformMatrix(startMatrix, startTranslation, startQuat, startScale);
	TransformConstraint::decomposeTransformMatrix(endMatrix, endTranslation, endQuat, endScale);

	// Interpolate translation
	//
	MVector translation = lerp(startTranslation, endTranslation, weight);
	MQuaternion quat = TransformConstraint::slerp(startQuat, endQuat, weight);
	MVector scale = lerp(startScale, endScale, weight);

	// Compose interpolated matrix
	//
	MMatrix translateMatrix = TransformConstraint::createPositionMatrix(translation);
	MMatrix rotateMatrix = quat.asMatrix();
	MMatrix scaleMatrix = TransformConstraint::createScaleMatrix(scale);

	return scaleMatrix * rotateMatrix * translateMatrix;

};


MMatrix TransformConstraint::blendMatrices(const MMatrix& restMatrix, const MMatrixArray& matrices, const MFloatArray& weights)
/**
Interpolates the supplied matrices using the weight array as a blend aplha.
The rest matrix is used just in case the weights don't equal 1.

@param restMatrix: The default matrix to blend from in case the weights don't equal 1.
@param matrices: The matrix array to blend.
@param weights: The float array containing the weighted averages.
@return: The interpolated matrix.
*/
{

	// Check the number of matrices
	//
	unsigned int numMatrices = matrices.length();

	switch (numMatrices)
	{

		case 0:
		{

			// Reuse rest matrix
			//
			return MMatrix(restMatrix);

		}
		break;

		case 1:
		{

			// Check weight sum before committing to either matrix
			//
			float weightSum = TransformConstraint::sum(weights);

			if (weightSum == 1.0f)
			{

				return MMatrix(matrices[0]);

			}
			else if (weightSum == 0.f)
			{

				return MMatrix(restMatrix);

			}
			else
			{

				return TransformConstraint::blendMatrices(restMatrix, matrices[0], weights[0]);

			}

		}
		break;

		default:
		{

			// Get start matrix
			//
			MFloatArray normalizedWeights = TransformConstraint::clamp(weights);
			float weightSum = TransformConstraint::sum(normalizedWeights);

			MMatrix matrix = MMatrix(matrices[0]);

			if (weightSum < 1.0f)
			{

				matrix = MMatrix(restMatrix);

			}

			// Get start transform components
			//
			unsigned int numMatrices = matrices.length();

			for (unsigned int i = 0; i < numMatrices; i++)
			{

				matrix = TransformConstraint::blendMatrices(matrix, matrices[i], normalizedWeights[i]);

			}

			return matrix;

		}
		break;

	}

	return MMatrix::identity;

};


double TransformConstraint::dot(const MQuaternion& quat, const MQuaternion& otherQuat)
/**
Returns the dot product of two quaternions.

@param quat: Quaternion.
@param: otherQuat: Other quaternion.
@return: Dot length.
*/
{

	return (quat.x * otherQuat.x) + (quat.y * otherQuat.y) + (quat.z * otherQuat.z) + (quat.w * otherQuat.w);

};


MQuaternion TransformConstraint::slerp(const MQuaternion& startQuat, const MQuaternion& endQuat, const float weight)
/**
Spherical interpolates two quaternions.
See the following for details: https://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/slerp/index.htm

@param startQuat: Start Quaternion.
@param endQuat: End Quaternion.
@param weight: The amount to interpolate.
@return: The interpolated quaternion.
*/
{

	MQuaternion q1 = MQuaternion(startQuat);
	MQuaternion q2 = MQuaternion(endQuat);

	double dot = TransformConstraint::dot(q1, q2);

	if (dot < 0.0)
	{

		dot = TransformConstraint::dot(q1, q2.negateIt());

	}

	double theta = acos(dot);
	double sinTheta = sin(theta);

	double w1, w2;

	if (sinTheta > 1e-3)
	{

		w1 = sin((1.0 - weight) * theta) / sinTheta;
		w2 = sin(weight * theta) / sinTheta;

	}
	else
	{

		w1 = 1.0 - weight;
		w2 = weight;

	}

	q1.scaleIt(w1);
	q2.scaleIt(w2);

	return q1 + q2;

};


float TransformConstraint::sum(const MFloatArray& items)
/**
Calculates the sum of all the supplied items.

@param items: The items to add up.
@return: The total sum.
*/
{

	// Iterate through numbers
	//
	unsigned int numItems = items.length();
	float sum = 0.0;

	for (unsigned int i = 0; i < numItems; i++)
	{

		sum += items[i];

	}

	return sum;

};


MFloatArray TransformConstraint::clamp(const MFloatArray& items)
/**
Clamps the supplied items so they don't exceed 1.
Anything below that is left alone and compensated for using the rest matrix.

@param items: The float array containing the weighted averages.
@return: The newly clamped array of weights.
*/
{

	// Check if sum is greater than one
	//
	float sum = TransformConstraint::sum(items);

	if (sum < 1.0)
	{

		return MFloatArray(items);

	}

	// Iterate through numbers
	//
	float fraction = 1.0f / sum;

	unsigned int numItems = items.length();
	MFloatArray normalizedItems = MFloatArray(numItems);

	for (unsigned int i = 0; i < numItems; i++)
	{

		normalizedItems[i] = items[i] * fraction;

	}

	return normalizedItems;


};


const MObject TransformConstraint::targetAttribute() const
/**
Returns the target attribute for the constraint.
Default implementation returns MObject::kNullObj.

@return: MObject
*/
{


	return TransformConstraint::target;

};


const MObject TransformConstraint::weightAttribute() const
/**
Returns the weight attribute for the constraint.
Default implementation returns MObject::kNullObj.

@return: MObject
*/
{


	return TransformConstraint::targetWeight;

};


const MObject TransformConstraint::constraintRotateOrderAttribute() const
/**
Returns the rotate order attribute for the constraint.
Default implementation returns MObject::kNullObj.

@return: MObject
*/
{


	return TransformConstraint::constraintRotateOrder;

};


void* TransformConstraint::creator()
/**
This function is called by Maya when a new instance is requested.
See pluginMain.cpp for details.

@return: TransformConstraint
*/
{

	return new TransformConstraint();

};


MStatus TransformConstraint::initialize()
/**
This function is called by Maya after a plugin has been loaded.
Use this function to define any static attributes.

@return: MStatus
*/
{

	MStatus status;

	// Declare attribute function sets
	//
	MFnCompoundAttribute fnCompoundAttr;
	MFnNumericAttribute fnNumericAttr;
	MFnUnitAttribute fnUnitAttr;
	MFnEnumAttribute fnEnumAttr;
	MFnMatrixAttribute fnMatrixAttr;

	// Input attributes:
	// ".restTranslateX" attribute
	//
	TransformConstraint::restTranslateX = fnUnitAttr.create("restTranslateX", "rtx", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::restCategory));

	// ".restTranslateY" attribute
	//
	TransformConstraint::restTranslateY = fnUnitAttr.create("restTranslateY", "rty", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::restCategory));

	// ".restTranslateZ" attribute
	//
	TransformConstraint::restTranslateZ = fnUnitAttr.create("restTranslateZ", "rtz", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::restCategory));

	// ".restTranslate" attribute
	//
	TransformConstraint::restTranslate = fnNumericAttr.create("restTranslate", "rt", TransformConstraint::restTranslateX, TransformConstraint::restTranslateY, TransformConstraint::restTranslateZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::restCategory));

	// ".restRotateX" attribute
	//
	TransformConstraint::restRotateX = fnUnitAttr.create("restRotateX", "rrx", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::restCategory));

	// ".restRotateY" attribute
	//
	TransformConstraint::restRotateY = fnUnitAttr.create("restRotateY", "rry", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::restCategory));

	// ".restRotateZ" attribute
	//
	TransformConstraint::restRotateZ = fnUnitAttr.create("restRotateZ", "rrz", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::restCategory));

	// ".restRotate" attribute
	//
	TransformConstraint::restRotate = fnNumericAttr.create("restRotate", "rr", TransformConstraint::restRotateX, TransformConstraint::restRotateY, TransformConstraint::restRotateZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::restCategory));

	// ".restScaleX" attribute
	//
	TransformConstraint::restScaleX = fnNumericAttr.create("restScaleX", "rsx", MFnNumericData::kDouble, 1.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::restCategory));

	// ".restScaleY" attribute
	//
	TransformConstraint::restScaleY = fnNumericAttr.create("restScaleY", "rsy", MFnNumericData::kDouble, 1.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::restCategory));

	// ".restScaleZ" attribute
	//
	TransformConstraint::restScaleZ = fnNumericAttr.create("restScaleZ", "rsz", MFnNumericData::kDouble, 1.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::restCategory));

	// ".restScale" attribute
	//
	TransformConstraint::restScale = fnNumericAttr.create("restScale", "rs", TransformConstraint::restScaleX, TransformConstraint::restScaleY, TransformConstraint::restScaleZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::restCategory));

	// ".constraintRotateOrder" attribute
	//
	TransformConstraint::constraintRotateOrder = fnEnumAttr.create("constraintRotateOrder", "cro", short(0), &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnEnumAttr.addField("xyz", 0));
	CHECK_MSTATUS(fnEnumAttr.addField("yzx", 1));
	CHECK_MSTATUS(fnEnumAttr.addField("zxy", 2));
	CHECK_MSTATUS(fnEnumAttr.addField("xzy", 3));
	CHECK_MSTATUS(fnEnumAttr.addField("yxz", 4));
	CHECK_MSTATUS(fnEnumAttr.addField("zyx", 5));
	CHECK_MSTATUS(fnEnumAttr.addToCategory(TransformConstraint::inputCategory));

	// ".constraintJointOrientX" attribute
	//
	TransformConstraint::constraintJointOrientX = fnUnitAttr.create("constraintJointOrientX", "cjox", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));

	// ".constraintJointOrientY" attribute
	//
	TransformConstraint::constraintJointOrientY = fnUnitAttr.create("constraintJointOrientY", "cjoy", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));

	// ".constraintJointOrientZ" attribute
	//
	TransformConstraint::constraintJointOrientZ = fnUnitAttr.create("constraintJointOrientZ", "cjoz", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));

	// ".constraintJointOrient" attribute
	//
	TransformConstraint::constraintJointOrient = fnNumericAttr.create("constraintJointOrient", "cjo", TransformConstraint::constraintJointOrientX, TransformConstraint::constraintJointOrientY, TransformConstraint::constraintJointOrientZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::inputCategory));

	// ".constraintParentInverseMatrix" attribute
	//
	TransformConstraint::constraintParentInverseMatrix = fnMatrixAttr.create("constraintParentInverseMatrix", "cpim", MFnMatrixAttribute::kDouble, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnMatrixAttr.addToCategory(TransformConstraint::inputCategory));

	// Target Attributes:
	// Define ".targetWeight" attribute
	//
	TransformConstraint::targetWeight = fnNumericAttr.create("targetWeight", "tw", MFnNumericData::kFloat, 0.0f, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.setMin(0.0));
	CHECK_MSTATUS(fnNumericAttr.setMax(1.0));

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetTranslateX" attribute
	//
	TransformConstraint::targetTranslateX = fnUnitAttr.create("targetTranslateX", "ttx", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetTranslateY" attribute
	//
	TransformConstraint::targetTranslateY = fnUnitAttr.create("targetTranslateY", "tty", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetTranslateZ" attribute
	//
	TransformConstraint::targetTranslateZ = fnUnitAttr.create("targetTranslateZ", "ttz", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetTranslate" attribute
	//
	TransformConstraint::targetTranslate = fnNumericAttr.create("targetTranslate", "tt", TransformConstraint::targetTranslateX, TransformConstraint::targetTranslateY, TransformConstraint::targetTranslateZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetOffsetTranslateX" attribute
	//
	TransformConstraint::targetOffsetTranslateX = fnUnitAttr.create("targetOffsetTranslateX", "totx", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetOffsetTranslateY" attribute
	//
	TransformConstraint::targetOffsetTranslateY = fnUnitAttr.create("targetOffsetTranslateY", "toty", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetOffsetTranslateZ" attribute
	//
	TransformConstraint::targetOffsetTranslateZ = fnUnitAttr.create("targetOffsetTranslateZ", "totz", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetOffsetTranslate" attribute
	//
	TransformConstraint::targetOffsetTranslate = fnNumericAttr.create("targetOffsetTranslate", "tot", TransformConstraint::targetOffsetTranslateX, TransformConstraint::targetOffsetTranslateY, TransformConstraint::targetOffsetTranslateZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetJointOrientX" attribute
	//
	TransformConstraint::targetJointOrientX = fnUnitAttr.create("targetJointOrientX", "tjox", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetJointOrientY" attribute
	//
	TransformConstraint::targetJointOrientY = fnUnitAttr.create("targetJointOrientY", "tjoy", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetJointOrientZ" attribute
	//
	TransformConstraint::targetJointOrientZ = fnUnitAttr.create("targetJointOrientZ", "tjoz", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetJointOrient" attribute
	//
	TransformConstraint::targetJointOrient = fnNumericAttr.create("targetJointOrient", "tjo", TransformConstraint::targetJointOrientX, TransformConstraint::targetJointOrientY, TransformConstraint::targetJointOrientZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetRotateX" attribute
	//
	TransformConstraint::targetRotateX = fnUnitAttr.create("targetRotateX", "trx", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetRotateY" attribute
	//
	TransformConstraint::targetRotateY = fnUnitAttr.create("targetRotateY", "try", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetRotateZ" attribute
	//
	TransformConstraint::targetRotateZ = fnUnitAttr.create("targetRotateZ", "trz", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetRotate" attribute
	//
	TransformConstraint::targetRotate = fnNumericAttr.create("targetRotate", "tr", TransformConstraint::targetRotateX, TransformConstraint::targetRotateY, TransformConstraint::targetRotateZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetOffsetRotateX" attribute
	//
	TransformConstraint::targetOffsetRotateX = fnUnitAttr.create("targetOffsetRotateX", "torx", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetOffsetRotateY" attribute
	//
	TransformConstraint::targetOffsetRotateY = fnUnitAttr.create("targetOffsetRotateY", "tory", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetOffsetRotateZ" attribute
	//
	TransformConstraint::targetOffsetRotateZ = fnUnitAttr.create("targetOffsetRotateZ", "torz", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetOffsetRotate" attribute
	//
	TransformConstraint::targetOffsetRotate = fnNumericAttr.create("targetOffsetRotate", "tor", TransformConstraint::targetOffsetRotateX, TransformConstraint::targetOffsetRotateY, TransformConstraint::targetOffsetRotateZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetRotateOrder" attribute
	//
	TransformConstraint::targetRotateOrder = fnEnumAttr.create("targetRotateOrder", "tro", short(0), &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnEnumAttr.addField("xyz", 0));
	CHECK_MSTATUS(fnEnumAttr.addField("yzx", 1));
	CHECK_MSTATUS(fnEnumAttr.addField("zxy", 2));
	CHECK_MSTATUS(fnEnumAttr.addField("xzy", 3));
	CHECK_MSTATUS(fnEnumAttr.addField("yxz", 4));
	CHECK_MSTATUS(fnEnumAttr.addField("zyx", 5));
	CHECK_MSTATUS(fnEnumAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnEnumAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetScaleX" attribute
	//
	TransformConstraint::targetScaleX = fnNumericAttr.create("targetScaleX", "tsx", MFnNumericData::kDouble, 1.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetScaleY" attribute
	//
	TransformConstraint::targetScaleY = fnNumericAttr.create("targetScaleY", "tsy", MFnNumericData::kDouble, 1.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetScaleZ" attribute
	//
	TransformConstraint::targetScaleZ = fnNumericAttr.create("targetScaleZ", "tsz", MFnNumericData::kDouble, 1.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetScale" attribute
	//
	TransformConstraint::targetScale = fnNumericAttr.create("targetScale", "ts", TransformConstraint::targetScaleX, TransformConstraint::targetScaleY, TransformConstraint::targetScaleZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetOffsetScaleX" attribute
	//
	TransformConstraint::targetOffsetScaleX = fnNumericAttr.create("targetOffsetScaleX", "tosx", MFnNumericData::kDouble, 1.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetOffsetScaleY" attribute
	//
	TransformConstraint::targetOffsetScaleY = fnNumericAttr.create("targetOffsetScaleY", "tosy", MFnNumericData::kDouble, 1.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetOffsetScaleZ" attribute
	//
	TransformConstraint::targetOffsetScaleZ = fnNumericAttr.create("targetOffsetScaleZ", "tosz", MFnNumericData::kDouble, 1.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetOffsetScale" attribute
	//
	TransformConstraint::targetOffsetScale = fnNumericAttr.create("targetOffsetScale", "tos", TransformConstraint::targetOffsetScaleX, TransformConstraint::targetOffsetScaleY, TransformConstraint::targetOffsetScaleZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetRotatePivotX" attribute
	//
	TransformConstraint::targetRotatePivotX = fnUnitAttr.create("targetRotatePivotX", "trpx", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetRotatePivotY" attribute
	//
	TransformConstraint::targetRotatePivotY = fnUnitAttr.create("targetRotatePivotY", "trpy", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetRotatePivotZ" attribute
	//
	TransformConstraint::targetRotatePivotZ = fnUnitAttr.create("targetRotatePivotZ", "trpz", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetRotatePivot" attribute
	//
	TransformConstraint::targetRotatePivot = fnNumericAttr.create("targetRotatePivot", "trp", TransformConstraint::targetRotatePivotX, TransformConstraint::targetRotatePivotY, TransformConstraint::targetRotatePivotZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetRotatePivotTranslateX" attribute
	//
	TransformConstraint::targetRotatePivotTranslateX = fnUnitAttr.create("targetRotatePivotTranslateX", "trptx", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetRotatePivotTranslateY" attribute
	//
	TransformConstraint::targetRotatePivotTranslateY = fnUnitAttr.create("targetRotatePivotTranslateY", "trpty", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetRotatePivotTranslateZ" attribute
	//
	TransformConstraint::targetRotatePivotTranslateZ = fnUnitAttr.create("targetRotatePivotTranslateZ", "trptz", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetRotatePivotTranslate" attribute
	//
	TransformConstraint::targetRotatePivotTranslate = fnNumericAttr.create("targetRotatePivotTranslate", "trpt", TransformConstraint::targetRotatePivotTranslateX, TransformConstraint::targetRotatePivotTranslateY, TransformConstraint::targetRotatePivotTranslateZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetScalePivotX" attribute
	//
	TransformConstraint::targetScalePivotX = fnUnitAttr.create("targetScalePivotX", "tspx", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetScalePivotY" attribute
	//
	TransformConstraint::targetScalePivotY = fnUnitAttr.create("targetScalePivotY", "tspy", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetScalePivotZ" attribute
	//
	TransformConstraint::targetScalePivotZ = fnUnitAttr.create("targetScalePivotZ", "tspz", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetScalePivot" attribute
	//
	TransformConstraint::targetScalePivot = fnNumericAttr.create("targetScalePivot", "tsp", TransformConstraint::targetScalePivotX, TransformConstraint::targetScalePivotY, TransformConstraint::targetScalePivotZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetScaleTranslateX" attribute
	//
	TransformConstraint::targetScalePivotTranslateX = fnUnitAttr.create("targetScalePivotTranslateX", "tsptx", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetScaleTranslateY" attribute
	//
	TransformConstraint::targetScalePivotTranslateY = fnUnitAttr.create("targetScalePivotTranslateY", "tspty", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetScaleTranslateZ" attribute
	//
	TransformConstraint::targetScalePivotTranslateZ = fnUnitAttr.create("targetScalePivotTranslateZ", "tsptz", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetScaleTranslate" attribute
	//
	TransformConstraint::targetScalePivotTranslate = fnNumericAttr.create("targetScalePivotTranslate", "tspt", TransformConstraint::targetScalePivotTranslateX, TransformConstraint::targetScalePivotTranslateY, TransformConstraint::targetScalePivotTranslateZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetParentMatrix" attribute
	//
	TransformConstraint::targetParentMatrix = fnMatrixAttr.create("targetParentMatrix", "tpm", MFnMatrixAttribute::kDouble, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnMatrixAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnMatrixAttr.addToCategory(TransformConstraint::targetCategory));

	// Define ".target" attribute
	//
	TransformConstraint::target = fnCompoundAttr.create("target", "tg",&status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnCompoundAttr.addChild(TransformConstraint::targetWeight));
	CHECK_MSTATUS(fnCompoundAttr.addChild(TransformConstraint::targetTranslate));
	CHECK_MSTATUS(fnCompoundAttr.addChild(TransformConstraint::targetOffsetTranslate));
	CHECK_MSTATUS(fnCompoundAttr.addChild(TransformConstraint::targetJointOrient));
	CHECK_MSTATUS(fnCompoundAttr.addChild(TransformConstraint::targetRotate));
	CHECK_MSTATUS(fnCompoundAttr.addChild(TransformConstraint::targetOffsetRotate));
	CHECK_MSTATUS(fnCompoundAttr.addChild(TransformConstraint::targetRotateOrder));
	CHECK_MSTATUS(fnCompoundAttr.addChild(TransformConstraint::targetScale));
	CHECK_MSTATUS(fnCompoundAttr.addChild(TransformConstraint::targetOffsetScale));
	CHECK_MSTATUS(fnCompoundAttr.addChild(TransformConstraint::targetRotatePivot));
	CHECK_MSTATUS(fnCompoundAttr.addChild(TransformConstraint::targetRotatePivotTranslate));
	CHECK_MSTATUS(fnCompoundAttr.addChild(TransformConstraint::targetScalePivot));
	CHECK_MSTATUS(fnCompoundAttr.addChild(TransformConstraint::targetScalePivotTranslate));
	CHECK_MSTATUS(fnCompoundAttr.addChild(TransformConstraint::targetParentMatrix));
	CHECK_MSTATUS(fnCompoundAttr.setArray(true));
	CHECK_MSTATUS(fnCompoundAttr.addToCategory(TransformConstraint::inputCategory));
	CHECK_MSTATUS(fnCompoundAttr.addToCategory(TransformConstraint::targetCategory));

	// Output attributes:
	// Define ".constraintTranslateX" attribute
	//
	TransformConstraint::constraintTranslateX = fnUnitAttr.create("constraintTranslateX", "ctx", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.setWritable(false));
	CHECK_MSTATUS(fnUnitAttr.setStorable(false));

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::outputCategory));

	// ".constraintTranslateY" attribute
	//
	TransformConstraint::constraintTranslateY = fnUnitAttr.create("constraintTranslateY", "cty", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.setWritable(false));
	CHECK_MSTATUS(fnUnitAttr.setStorable(false));

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::outputCategory));

	// ".constraintTranslateZ" attribute
	//
	TransformConstraint::constraintTranslateZ = fnUnitAttr.create("constraintTranslateZ", "ctz", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.setWritable(false));
	CHECK_MSTATUS(fnUnitAttr.setStorable(false));

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::outputCategory));

	// ".constraintTranslate" attribute
	//
	TransformConstraint::constraintTranslate = fnNumericAttr.create("constraintTranslate", "ct", TransformConstraint::constraintTranslateX, TransformConstraint::constraintTranslateY, TransformConstraint::constraintTranslateZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.setWritable(false));
	CHECK_MSTATUS(fnNumericAttr.setStorable(false));

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::outputCategory));

	// ".constraintRotateX" attribute
	//
	TransformConstraint::constraintRotateX = fnUnitAttr.create("constraintRotateX", "crx", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.setWritable(false));
	CHECK_MSTATUS(fnUnitAttr.setStorable(false));

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::outputCategory));

	// ".constraintRotateY" attribute
	//
	TransformConstraint::constraintRotateY = fnUnitAttr.create("constraintRotateY", "cry", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.setWritable(false));
	CHECK_MSTATUS(fnUnitAttr.setStorable(false));

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::outputCategory));

	// ".constraintRotateZ" attribute
	//
	TransformConstraint::constraintRotateZ = fnUnitAttr.create("constraintRotateZ", "crz", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.setWritable(false));
	CHECK_MSTATUS(fnUnitAttr.setStorable(false));

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::outputCategory));

	// ".constraintRotate" attribute
	//
	TransformConstraint::constraintRotate = fnNumericAttr.create("constraintRotate", "cr", TransformConstraint::constraintRotateX, TransformConstraint::constraintRotateY, TransformConstraint::constraintRotateZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.setWritable(false));
	CHECK_MSTATUS(fnNumericAttr.setStorable(false));

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::outputCategory));

	// ".constraintScaleX" attribute
	//
	TransformConstraint::constraintScaleX = fnNumericAttr.create("constraintScaleX", "csx", MFnNumericData::kDouble, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.setWritable(false));
	CHECK_MSTATUS(fnNumericAttr.setStorable(false));

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::outputCategory));

	// ".constraintScaleY" attribute
	//
	TransformConstraint::constraintScaleY = fnNumericAttr.create("constraintScaleY", "csy", MFnNumericData::kDouble, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.setWritable(false));
	CHECK_MSTATUS(fnNumericAttr.setStorable(false));

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::outputCategory));

	// ".constraintScaleZ" attribute
	//
	TransformConstraint::constraintScaleZ = fnNumericAttr.create("constraintScaleZ", "csz", MFnNumericData::kDouble, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.setWritable(false));
	CHECK_MSTATUS(fnNumericAttr.setStorable(false));

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::outputCategory));

	// ".constraintScale" attribute
	//
	TransformConstraint::constraintScale = fnNumericAttr.create("constraintScale", "cs", TransformConstraint::constraintScaleX, TransformConstraint::constraintScaleY, TransformConstraint::constraintScaleZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.setWritable(false));
	CHECK_MSTATUS(fnNumericAttr.setStorable(false));

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::outputCategory));

	// ".constraintMatrix" attribute
	//
	TransformConstraint::constraintMatrix = fnMatrixAttr.create("constraintMatrix", "cm", MFnMatrixAttribute::kDouble, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnMatrixAttr.setWritable(false));
	CHECK_MSTATUS(fnMatrixAttr.setStorable(false));

	CHECK_MSTATUS(fnMatrixAttr.addToCategory(TransformConstraint::outputCategory));

	// ".constraintInverseMatrix" attribute
	//
	TransformConstraint::constraintInverseMatrix = fnMatrixAttr.create("constraintInverseMatrix", "cim", MFnMatrixAttribute::kDouble, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnMatrixAttr.setWritable(false));
	CHECK_MSTATUS(fnMatrixAttr.setStorable(false));

	CHECK_MSTATUS(fnMatrixAttr.addToCategory(TransformConstraint::outputCategory));

	// ".constraintWorldMatrix" attribute
	//
	TransformConstraint::constraintWorldMatrix = fnMatrixAttr.create("constraintWorldMatrix", "cwm", MFnMatrixAttribute::kDouble, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnMatrixAttr.setWritable(false));
	CHECK_MSTATUS(fnMatrixAttr.setStorable(false));

	CHECK_MSTATUS(fnMatrixAttr.addToCategory(TransformConstraint::outputCategory));

	// ".constraintWorldInverseMatrix" attribute
	//
	TransformConstraint::constraintWorldInverseMatrix = fnMatrixAttr.create("constraintWorldInverseMatrix", "cwim", MFnMatrixAttribute::kDouble, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnMatrixAttr.setWritable(false));
	CHECK_MSTATUS(fnMatrixAttr.setStorable(false));

	CHECK_MSTATUS(fnMatrixAttr.addToCategory(TransformConstraint::outputCategory));

	// Add attributes
	//
	CHECK_MSTATUS(TransformConstraint::addAttribute(TransformConstraint::restTranslate));
	CHECK_MSTATUS(TransformConstraint::addAttribute(TransformConstraint::restRotate));
	CHECK_MSTATUS(TransformConstraint::addAttribute(TransformConstraint::restScale));
	CHECK_MSTATUS(TransformConstraint::addAttribute(TransformConstraint::target));
	CHECK_MSTATUS(TransformConstraint::addAttribute(TransformConstraint::constraintTranslate));
	CHECK_MSTATUS(TransformConstraint::addAttribute(TransformConstraint::constraintRotate));
	CHECK_MSTATUS(TransformConstraint::addAttribute(TransformConstraint::constraintRotateOrder));
	CHECK_MSTATUS(TransformConstraint::addAttribute(TransformConstraint::constraintJointOrient));
	CHECK_MSTATUS(TransformConstraint::addAttribute(TransformConstraint::constraintScale));
	CHECK_MSTATUS(TransformConstraint::addAttribute(TransformConstraint::constraintMatrix));
	CHECK_MSTATUS(TransformConstraint::addAttribute(TransformConstraint::constraintInverseMatrix));
	CHECK_MSTATUS(TransformConstraint::addAttribute(TransformConstraint::constraintWorldMatrix));
	CHECK_MSTATUS(TransformConstraint::addAttribute(TransformConstraint::constraintWorldInverseMatrix));
	CHECK_MSTATUS(TransformConstraint::addAttribute(TransformConstraint::constraintParentInverseMatrix));

	// Define target attribute relationships
	//
	status = fnCompoundAttr.setObject(TransformConstraint::target);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	unsigned int numChildren = fnCompoundAttr.numChildren();

	for (unsigned int i = 0; i < numChildren; i++) {

		MObject child = fnCompoundAttr.child(i);

		CHECK_MSTATUS(TransformConstraint::attributeAffects(child, TransformConstraint::constraintTranslate));
		CHECK_MSTATUS(TransformConstraint::attributeAffects(child, TransformConstraint::constraintTranslateX));
		CHECK_MSTATUS(TransformConstraint::attributeAffects(child, TransformConstraint::constraintTranslateY));
		CHECK_MSTATUS(TransformConstraint::attributeAffects(child, TransformConstraint::constraintTranslateZ));

		CHECK_MSTATUS(TransformConstraint::attributeAffects(child, TransformConstraint::constraintRotate));
		CHECK_MSTATUS(TransformConstraint::attributeAffects(child, TransformConstraint::constraintRotateX));
		CHECK_MSTATUS(TransformConstraint::attributeAffects(child, TransformConstraint::constraintRotateY));
		CHECK_MSTATUS(TransformConstraint::attributeAffects(child, TransformConstraint::constraintRotateZ));

		CHECK_MSTATUS(TransformConstraint::attributeAffects(child, TransformConstraint::constraintScale));
		CHECK_MSTATUS(TransformConstraint::attributeAffects(child, TransformConstraint::constraintScaleX));
		CHECK_MSTATUS(TransformConstraint::attributeAffects(child, TransformConstraint::constraintScaleY));
		CHECK_MSTATUS(TransformConstraint::attributeAffects(child, TransformConstraint::constraintScaleZ));

		CHECK_MSTATUS(TransformConstraint::attributeAffects(child, TransformConstraint::constraintMatrix));
		CHECK_MSTATUS(TransformConstraint::attributeAffects(child, TransformConstraint::constraintInverseMatrix));
		CHECK_MSTATUS(TransformConstraint::attributeAffects(child, TransformConstraint::constraintWorldMatrix));
		CHECK_MSTATUS(TransformConstraint::attributeAffects(child, TransformConstraint::constraintWorldInverseMatrix));

	}

	// Define rest attribute relationships
	//
	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::restTranslate, TransformConstraint::constraintTranslate));
	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::restTranslateX, TransformConstraint::constraintTranslateX));
	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::restTranslateY, TransformConstraint::constraintTranslateY));
	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::restTranslateZ, TransformConstraint::constraintTranslateZ));

	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::restRotate, TransformConstraint::constraintRotate));
	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::restRotateX, TransformConstraint::constraintRotateX));
	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::restRotateY, TransformConstraint::constraintRotateY));
	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::restRotateZ, TransformConstraint::constraintRotateZ));

	// Define constraint attribute relationships
	//
	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::constraintRotateOrder, TransformConstraint::constraintRotate));

	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::constraintJointOrient, TransformConstraint::constraintRotate));
	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::constraintJointOrientX, TransformConstraint::constraintRotateX));
	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::constraintJointOrientY, TransformConstraint::constraintRotateY));
	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::constraintJointOrientZ, TransformConstraint::constraintRotateZ));

	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::constraintParentInverseMatrix, TransformConstraint::constraintTranslate));
	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::constraintParentInverseMatrix, TransformConstraint::constraintTranslateX));
	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::constraintParentInverseMatrix, TransformConstraint::constraintTranslateY));
	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::constraintParentInverseMatrix, TransformConstraint::constraintTranslateZ));
	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::constraintParentInverseMatrix, TransformConstraint::constraintRotate));
	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::constraintParentInverseMatrix, TransformConstraint::constraintRotateX));
	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::constraintParentInverseMatrix, TransformConstraint::constraintRotateY));
	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::constraintParentInverseMatrix, TransformConstraint::constraintRotateZ));
	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::constraintParentInverseMatrix, TransformConstraint::constraintScale));
	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::constraintParentInverseMatrix, TransformConstraint::constraintScaleX));
	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::constraintParentInverseMatrix, TransformConstraint::constraintScaleY));
	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::constraintParentInverseMatrix, TransformConstraint::constraintScaleZ));
	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::constraintParentInverseMatrix, TransformConstraint::constraintMatrix));
	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::constraintParentInverseMatrix, TransformConstraint::constraintInverseMatrix));
	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::constraintParentInverseMatrix, TransformConstraint::constraintWorldMatrix));
	CHECK_MSTATUS(TransformConstraint::attributeAffects(TransformConstraint::constraintParentInverseMatrix, TransformConstraint::constraintWorldInverseMatrix));

	return MS::kSuccess;

};
