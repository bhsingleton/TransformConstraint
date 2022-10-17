#ifndef _TransformConstraintNode
#define _TransformConstraintNode
//
// File: TransformConstraintNode.h
//
// Dependency Graph Node: TransformConstraint
//
// Author: Ben Singleton
//

#include <maya/MPxConstraint.h>

#include <maya/MGlobal.h>
#include <maya/MTypeId.h>
#include <maya/MObject.h>
#include <maya/MObjectArray.h>
#include <maya/MObjectHandle.h>
#include <maya/MDagPath.h>
#include <maya/MUuid.h>
#include <maya/MPlug.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MArrayDataHandle.h>
#include <maya/MDistance.h>
#include <maya/MAngle.h>
#include <maya/MVector.h>
#include <maya/MEulerRotation.h>
#include <maya/MQuaternion.h>
#include <maya/MFloatArray.h>
#include <maya/MMatrix.h>
#include <maya/MMatrixArray.h>
#include <maya/MFnCompoundAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnNumericData.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnMessageAttribute.h>
#include <maya/MFnData.h>
#include <maya/MFnMatrixData.h>
#include <maya/MDagMessage.h>
#include <maya/MSelectionList.h>
#include <maya/MItDependencyNodes.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MFnReference.h>
#include <maya/MDagModifier.h>
#include <maya/MCallbackIdArray.h>

#include "map"


class TransformConstraint : public MPxConstraint
{

public:

							TransformConstraint();
	virtual					~TransformConstraint(); 

	virtual MStatus			compute(const MPlug& plug, MDataBlock& data);

	static  void*			creator();
	static  MStatus			initialize();

	const	MObject			targetAttribute() const override;
	const	MObject			weightAttribute() const override;
	const	MObject			constraintRotateOrderAttribute() const override;
	
	static	MStatus			decomposeTransformMatrix(const MMatrix& matrix, MVector& position, MQuaternion& rotation, MVector& scale);
	static	MVector			getTranslationPart(const MMatrix& matrix);
	static	MQuaternion		getRotationPart(const MMatrix& matrix);
	static	MVector			getScalePart(const MMatrix& matrix);

	static	MMatrix			createPositionMatrix(const MVector& position);
	static	MMatrix			createPositionMatrix(const MMatrix& matrix);

	static	MMatrix			createRotationMatrix(const MVector& rotation, const int rotateOrder);
	static	MMatrix			createRotationMatrix(const MMatrix& matrix);

	static	MMatrix			createScaleMatrix(const MVector& scale);
	static	MMatrix			createScaleMatrix(const MMatrix& matrix);

	static	double			dot(const MQuaternion& startQuat, const MQuaternion& endQuat);
	static	MQuaternion		slerp(const MQuaternion& startQuat, const MQuaternion& endQuat, const float weight);

	static	MMatrix			blendMatrices(const MMatrix& restMatrix, const MMatrixArray& matrices, const MFloatArray& weights);
	static	MMatrix			blendMatrices(const MMatrix& startMatrix, const MMatrix& endMatrix, const float weight);

	static	float			sum(const MFloatArray& items);
	static	MFloatArray		clamp(const MFloatArray& items);

public:

	static	MObject		restTranslate;
	static	MObject		restTranslateX;
	static	MObject		restTranslateY;
	static	MObject		restTranslateZ;
	static	MObject		restRotate;
	static	MObject		restRotateX;
	static	MObject		restRotateY;
	static	MObject		restRotateZ;
	static	MObject		restScale;
	static	MObject		restScaleX;
	static	MObject		restScaleY;
	static	MObject		restScaleZ;

	static	MObject		target;
	static	MObject		targetWeight;
	static	MObject		targetParentMatrix;
	static	MObject		targetTranslate;
	static	MObject		targetTranslateX;
	static	MObject		targetTranslateY;
	static	MObject		targetTranslateZ;
	static	MObject		targetOffsetTranslate;
	static	MObject		targetOffsetTranslateX;
	static	MObject		targetOffsetTranslateY;
	static	MObject		targetOffsetTranslateZ;
	static	MObject		targetJointOrient;
	static	MObject		targetJointOrientX;
	static	MObject		targetJointOrientY;
	static	MObject		targetJointOrientZ;
	static	MObject		targetRotate;
	static	MObject		targetRotateX;
	static	MObject		targetRotateY;
	static	MObject		targetRotateZ;
	static	MObject		targetOffsetRotate;
	static	MObject		targetOffsetRotateX;
	static	MObject		targetOffsetRotateY;
	static	MObject		targetOffsetRotateZ;
	static	MObject		targetRotateOrder;
	static	MObject		targetScale;
	static	MObject		targetScaleX;
	static	MObject		targetScaleY;
	static	MObject		targetScaleZ;
	static	MObject		targetOffsetScale;
	static	MObject		targetOffsetScaleX;
	static	MObject		targetOffsetScaleY;
	static	MObject		targetOffsetScaleZ;
	static	MObject		targetRotatePivot;
	static	MObject		targetRotatePivotX;
	static	MObject		targetRotatePivotY;
	static	MObject		targetRotatePivotZ;
	static	MObject		targetRotateTranslate;
	static	MObject		targetRotateTranslateX;
	static	MObject		targetRotateTranslateY;
	static	MObject		targetRotateTranslateZ;
	static	MObject		targetScalePivot;
	static	MObject		targetScalePivotX;
	static	MObject		targetScalePivotY;
	static	MObject		targetScalePivotZ;
	static	MObject		targetScaleTranslate;
	static	MObject		targetScaleTranslateX;
	static	MObject		targetScaleTranslateY;
	static	MObject		targetScaleTranslateZ;

	static	MObject		constraintTranslate;
	static	MObject		constraintTranslateX;
	static	MObject		constraintTranslateY;
	static	MObject		constraintTranslateZ;

	static	MObject		constraintRotate;
	static	MObject		constraintRotateX;
	static	MObject		constraintRotateY;
	static	MObject		constraintRotateZ;
	static	MObject		constraintRotateOrder;

	static	MObject		constraintScale;
	static	MObject		constraintScaleX;
	static	MObject		constraintScaleY;
	static	MObject		constraintScaleZ;

	static	MObject		constraintJointOrient;
	static	MObject		constraintJointOrientX;
	static	MObject		constraintJointOrientY;
	static	MObject		constraintJointOrientZ;
	static	MObject		constraintMatrix;
	static	MObject		constraintInverseMatrix;
	static	MObject		constraintWorldMatrix;
	static	MObject		constraintWorldInverseMatrix;
	static	MObject		constraintParentInverseMatrix;
	static	MObject		constraintObject;

public:

	static	MTypeId		id;
	static	MString		targetCategory;
	static	MString		outputCategory;

};

#endif
