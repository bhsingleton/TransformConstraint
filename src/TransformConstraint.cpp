//
// File: TransformConstraintNode.cpp
//
// Dependency Graph Node: TransformConstraint
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
MObject TransformConstraint::targetRotateTranslate;
MObject TransformConstraint::targetRotateTranslateX;
MObject TransformConstraint::targetRotateTranslateY;
MObject TransformConstraint::targetRotateTranslateZ;
MObject TransformConstraint::targetScalePivot;
MObject TransformConstraint::targetScalePivotX;
MObject TransformConstraint::targetScalePivotY;
MObject TransformConstraint::targetScalePivotZ;
MObject TransformConstraint::targetScaleTranslate;
MObject TransformConstraint::targetScaleTranslateX;
MObject TransformConstraint::targetScaleTranslateY;
MObject TransformConstraint::targetScaleTranslateZ;

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
MObject TransformConstraint::constraintObject;

MTypeId	TransformConstraint::id(0x00131804);
MString	TransformConstraint::targetCategory("Target");
MString	TransformConstraint::outputCategory("Output");

std::map<long, TransformConstraint*> TransformConstraint::instances = std::map<long, TransformConstraint*>();
MCallbackId	TransformConstraint::childAddedCallbackId;


TransformConstraint::TransformConstraint() {}


TransformConstraint::~TransformConstraint()
/**
Destructor.
*/
{

	this->instances.erase(this->hashCode());

};


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
		int constraintRotateOrder = constraintRotateOrderHandle.asInt();

		MVector restTranslate = restTranslateHandle.asVector();
		MEulerRotation restRotate = MEulerRotation(restRotateHandle.asVector(), MEulerRotation::RotationOrder(constraintRotateOrder));
		double3 &restScale = restScaleHandle.asDouble3();

		MMatrix restTranslateMatrix = TransformConstraint::createPositionMatrix(restTranslate);
		MMatrix restRotateMatrix = restRotate.asMatrix();
		MMatrix restScaleMatrix = TransformConstraint::createScaleMatrix(restScale);

		MMatrix constraintParentInverseMatrix = constraintParentInverseMatrixHandle.asMatrix();

		MMatrix restMatrix = restScaleMatrix * restRotateMatrix * restTranslateMatrix;
		MMatrix restWorldMatrix = restMatrix * constraintParentInverseMatrix.inverse();

		MEulerRotation constraintJointOrient = MEulerRotation(constraintJointOrientHandle.asVector(), MEulerRotation::RotationOrder(constraintRotateOrder));

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
		MDataHandle targetRotatePivotHandle, targetRotateTranslateHandle;
		MDataHandle targetScalePivotHandle, targetScaleTranslateHandle;

		int targetRotateOrder;
		MMatrix targetOffsetMatrix;
		MMatrix targetParentMatrix;
		MMatrix targetTranslateMatrix, targetOffsetTranslateMatrix;
		MMatrix targetJointOrientMatrix, targetRotateMatrix, targetOffsetRotateMatrix;
		MMatrix targetScaleMatrix, targetOffsetScaleMatrix;
		MMatrix targetRotatePivotMatrix, targetRotateTranslateMatrix;
		MMatrix	targetScalePivotMatrix, targetScaleTranslateMatrix;

		for (unsigned int i = 0; i < targetCount; i++)
		{

			// Jump to array element
			//
			status = targetArrayHandle.jumpToElement(i);
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
			targetRotateTranslateHandle = targetHandle.child(TransformConstraint::targetRotateTranslate);
			targetScalePivotHandle = targetHandle.child(TransformConstraint::targetScalePivot);
			targetScaleTranslateHandle = targetHandle.child(TransformConstraint::targetScaleTranslate);

			// Get weight value
			//
			targetWeights[i] = targetWeightHandle.asFloat();

			// Get rotate order
			//
			targetRotateOrder = targetRotateOrderHandle.asInt();

			// Compute target offset matrix
			//
			targetOffsetTranslateMatrix = TransformConstraint::createPositionMatrix(targetOffsetTranslateHandle.asVector());
			targetOffsetRotateMatrix = TransformConstraint::createRotationMatrix(targetOffsetRotateHandle.asVector(), targetRotateOrder);
			targetOffsetScaleMatrix = TransformConstraint::createScaleMatrix(targetOffsetScaleHandle.asDouble3());

			targetOffsetMatrix = targetOffsetScaleMatrix * targetOffsetRotateMatrix * targetOffsetTranslateMatrix;

			// Compute target matrices
			//
			targetParentMatrix = targetParentMatrixHandle.asMatrix();
			targetTranslateMatrix = TransformConstraint::createPositionMatrix(targetTranslateHandle.asVector());
			
			targetJointOrientMatrix = TransformConstraint::createRotationMatrix(targetJointOrientHandle.asVector(), targetRotateOrder);
			targetRotateMatrix = TransformConstraint::createRotationMatrix(targetRotateHandle.asVector(), targetRotateOrder);
			
			targetScaleMatrix = TransformConstraint::createScaleMatrix(targetScaleHandle.asDouble3());
			
			targetRotatePivotMatrix = TransformConstraint::createPositionMatrix(targetRotatePivotHandle.asVector());
			targetRotateTranslateMatrix = TransformConstraint::createPositionMatrix(targetRotateTranslateHandle.asVector());
			targetScalePivotMatrix = TransformConstraint::createPositionMatrix(targetScalePivotHandle.asVector());
			targetScaleTranslateMatrix = TransformConstraint::createPositionMatrix(targetScaleTranslateHandle.asVector());

			targetMatrices[i] = targetScalePivotMatrix * targetScaleMatrix * targetScaleTranslateMatrix * targetRotatePivotMatrix * targetRotateMatrix * targetJointOrientMatrix * targetRotateTranslateMatrix * targetTranslateMatrix;
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

		// Initialize transformation matrix
		//
		MTransformationMatrix transformationMatrix = MTransformationMatrix(MMatrix(constraintMatrix));

		status = transformationMatrix.reorderRotation(MTransformationMatrix::RotationOrder(constraintRotateOrder + 1));
		CHECK_MSTATUS_AND_RETURN_IT(status);

		// Compensate for any joint orientation
		//
		transformationMatrix.rotateBy(constraintJointOrient.inverse(), MSpace::kTransform, &status);
		CHECK_MSTATUS_AND_RETURN_IT(status);

		// Set translation constraint
		//
		MVector constraintTranslate = transformationMatrix.getTranslation(MSpace::kTransform, &status);
		CHECK_MSTATUS_AND_RETURN_IT(status);

		MDistance::Unit internalUnit = MDistance::internalUnit();

		constraintTranslateXHandle.setMDistance(MDistance(constraintTranslate.x, internalUnit));
		constraintTranslateYHandle.setMDistance(MDistance(constraintTranslate.y, internalUnit));
		constraintTranslateZHandle.setMDistance(MDistance(constraintTranslate.z, internalUnit));

		constraintTranslateXHandle.setClean();
		constraintTranslateYHandle.setClean();
		constraintTranslateZHandle.setClean();

		// Set rotation constraint
		//
		MEulerRotation constraintRotate = transformationMatrix.eulerRotation();

		constraintRotateXHandle.setMAngle(MAngle(constraintRotate.x, MAngle::kRadians));
		constraintRotateYHandle.setMAngle(MAngle(constraintRotate.y, MAngle::kRadians));
		constraintRotateZHandle.setMAngle(MAngle(constraintRotate.z, MAngle::kRadians));

		constraintRotateXHandle.setClean();
		constraintRotateYHandle.setClean();
		constraintRotateZHandle.setClean();

		// Set scale constraint
		//
		double constraintScale[3];
		status = transformationMatrix.getScale(constraintScale, MSpace::kTransform);

		constraintScaleXHandle.setDouble(constraintScale[0]);
		constraintScaleYHandle.setDouble(constraintScale[1]);
		constraintScaleZHandle.setDouble(constraintScale[2]);

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
	else;

	return MS::kUnknownParameter;

}


void onChildAdded(MDagPath& child, MDagPath& parent, void* clientData)
/**
Child added callback function.
This function will handle updating the constraint's parent handle.

@param child: The newly added child.
@param parent: The parent owner.
@param clientData: Pointer to any client data passed on creation.
@return: void
*/
{

	MStatus status;

	// Iterate through instances
	//
	std::map<long, TransformConstraint*>::iterator iter;

	long hashCode;
	TransformConstraint* constraint;

	for (iter = TransformConstraint::instances.begin(); iter != TransformConstraint::instances.end(); iter++)
	{

		// Check if child is derived from constraint
		//
		hashCode = iter->first;
		constraint = iter->second;

		if (constraint->constraintHandle.object() == child.node())
		{

			status = constraint->updateConstraintParentInverseMatrix();
			CHECK_MSTATUS(status);

		}

	}

};


MStatus TransformConstraint::legalConnection(const MPlug& plug, const MPlug& otherPlug, bool asSrc, bool& isLegal)
/**
This method allows you to check for legal connections being made to attributes of this node.
You should return kUnknownParameter to specify that maya should handle this connection if you are unable to determine if it is legal.

@param plug: Attribute on this node.
@param otherPlug: Attribute on other node.
@param asSrc: Is this plug a source of the connection.
@param isLegal: Set this to true if the connection is legal otherwise false.
@return: MStatus
*/
{

	// Check the plug attribute
	//
	MObject attribute = plug.attribute();
	
	if (attribute == TransformConstraint::constraintObject && !asSrc)
	{

		// Verify other node is a dag node
		//
		MObject otherNode = otherPlug.node();
		MObject otherAttribute = otherPlug.attribute();

		if (otherNode.hasFn(MFn::kDagNode) && otherAttribute.hasFn(MFn::kMessageAttribute))
		{

			isLegal = true;

		}
		else
		{

			isLegal = false;

		}

		return MS::kSuccess;

	}

	return MS::kUnknownParameter;

};


MStatus TransformConstraint::connectionMade(const MPlug& plug, const MPlug& otherPlug, bool asSrc)
/**
This method gets called when connections are made to attributes of this node.
You should return kUnknownParameter to specify that maya should handle this connection or if you want maya to process the connection as well.

@param plug: Attribute on this node.
@param otherPlug: Attribute on other node.
@param asSrc: Is this plug a source of the connection.
@return: MStatus
*/
{

	MStatus status;

	// Check the plug attribute
	//
	MObject attribute = plug.attribute();

	if (attribute == TransformConstraint::constraintObject && !asSrc)
	{

		this->constraintHandle = MObjectHandle(otherPlug.node());

	}

	return MS::kUnknownParameter;

};


MStatus TransformConstraint::connectionBroken(const MPlug& plug, const MPlug& otherPlug, bool asSrc)
/**
This method gets called when connections are broken with attributes of this node.
You should return kUnknownParameter to specify that maya should handle this connection or if you want maya to process the connection as well.

@param plug: Attribute on this node.
@param otherPlug: Attribute on other node.
@param asSrc: Is this plug a source of the connection.
@return: MStatus
*/
{

	MStatus status;

	// Check the plug attribute
	//
	MObject attribute = plug.attribute();

	if (attribute == TransformConstraint::constraintObject && !asSrc)
	{

		this->constraintHandle = MObjectHandle();

	}

	return MS::kUnknownParameter;

};


MStatus TransformConstraint::connectPlugs(MPlug& source, MPlug& destination)
/**
Connects the two supplied plugs.

@param source: The source plug.
@param destination: The destination plug.
@return: Return status;
*/
{

	MStatus status;

	// Check if plugs are valid
	//
	if (source.isNull() || destination.isNull())
	{

		return MS::kFailure;

	}

	// Execute dag modifier
	//
	MDagModifier dagModifier;

	status = dagModifier.connect(source, destination);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	return dagModifier.doIt();

};


MStatus TransformConstraint::disconnectPlugs(MPlug& source, MPlug& destination)
/**
Disconnects the two supplied plugs.

@param source: The source plug.
@param destination: The destination plug.
@return: Return status;
*/
{

	MStatus status;

	// Check if plugs are valid
	//
	if (source.isNull() || destination.isNull())
	{

		return MS::kFailure;

	}

	// Execute dag modifier
	//
	MDagModifier dagModifier;

	status = dagModifier.disconnect(source, destination);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	return dagModifier.doIt();

};


MStatus TransformConstraint::breakConnections(MPlug& plug, bool source, bool destination)
/**
Breaks the specified connections from the supplied plug.

@param plug: The plug the break connections from.
@param source: Specifies if the source plug should be broken.
@param destination: Specifies if the destination plugs should be broken.
@return: Return status;
*/
{

	MStatus status;

	// Check if source connections should be broken
	//
	if (source)
	{

		MPlug otherPlug = plug.source();

		if (!otherPlug.isNull())
		{

			status = TransformConstraint::disconnectPlugs(otherPlug, plug);
			CHECK_MSTATUS_AND_RETURN_IT(status);

		}

	}

	// Check if destination plugs should be broken
	//
	if (destination)
	{

		MPlugArray otherPlugs;

		bool isConnected = plug.destinations(otherPlugs, &status);
		CHECK_MSTATUS_AND_RETURN_IT(status);

		unsigned int numOtherPlugs = otherPlugs.length();

		for (unsigned int i = 0; i < numOtherPlugs; i++)
		{

			status = TransformConstraint::disconnectPlugs(plug, otherPlugs[i]);
			CHECK_MSTATUS_AND_RETURN_IT(status);

		}

	}

	return status;

};


MStatus TransformConstraint::updateConstraintParentInverseMatrix()
/**
Updates the plug connected to the constraintParentInverseMatrix plug.
This function should be called by the childAdded callback whenever the constraint object's parent is changed.
This will ensure the correct worldMatrix plug is used.

@return: Return status.

*/
{

	MStatus status;

	// Check if constraint handle is still alive
	//
	if (!this->constraintHandle.isAlive())
	{

		return MS::kFailure;

	}

	MObject constraintObject = this->constraintHandle.object();

	// Break connections to plug
	//
	MPlug constraintParentInverseMatrixPlug = MPlug(this->thisMObject(), TransformConstraint::constraintParentInverseMatrix);

	status = TransformConstraint::breakConnections(constraintParentInverseMatrixPlug, true, false);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	// Check if constraint object has a parent
	//
	MObject parent = TransformConstraint::getParentOf(constraintObject);

	if (!parent.isNull())
	{

		// Initialize function set from dag path
		//
		MDagPath dagPath = TransformConstraint::getAPathTo(parent);

		MFnDagNode fnDagNode(dagPath, &status);
		CHECK_MSTATUS_AND_RETURN_IT(status);

		// Connect parent's worldInverseMatrix plug to constraintParentInverseMatrix plug
		// This will bypass any cyclical dependencies from the offsetParentMatrix plug!
		//
		MPlug worldInverseMatrixPlug = fnDagNode.findPlug("worldInverseMatrix", false, &status);
		CHECK_MSTATUS_AND_RETURN_IT(status);

		status = worldInverseMatrixPlug.selectAncestorLogicalIndex(dagPath.instanceNumber());
		CHECK_MSTATUS_AND_RETURN_IT(status);

		status = TransformConstraint::connectPlugs(worldInverseMatrixPlug, constraintParentInverseMatrixPlug);
		CHECK_MSTATUS_AND_RETURN_IT(status);

	}
	else
	{

		// Reset plug matrix
		//
		MObject matrixData = TransformConstraint::createMatrixData(MMatrix::identity);

		status = constraintParentInverseMatrixPlug.setMObject(matrixData);
		CHECK_MSTATUS_AND_RETURN_IT(status);

	}

	return status;

};


MObject TransformConstraint::createMatrixData(MMatrix matrix)
/**
Converts a matrix into a data object compatible with plugs.

@param matrix: The matrix to convert.
@return: A matrix data object.
*/
{

	MFnMatrixData fnMatrixData;

	MObject matrixData = fnMatrixData.create();
	fnMatrixData.set(matrix);

	return matrixData;

};


MMatrix TransformConstraint::getMatrixData(MObject& matrixData)
/**
Converts a data object into a matrix compatible with plugs.

@param matrixData: The data object to convert.
@return: A matrix.
*/
{

	MFnMatrixData fnMatrixData(matrixData);
	return fnMatrixData.matrix();

};


MDagPath TransformConstraint::getAPathTo(MObject& dependNode)
/**
Returns the dag path for the supplied dag node.

@param dependNode: The node to get a path to.
@return: The dag path.
*/
{

	// Verify this is a dag node
	//
	if (!dependNode.hasFn(MFn::kDagNode))
	{

		return MDagPath();

	}

	// Get a path to dag node
	//
	MDagPath dagPath;
	MDagPath::getAPathTo(dependNode, dagPath);

	return dagPath;

}


MObject TransformConstraint::getParentOf(MObject& dependNode)
/**
Returns the parent for the supplied dag node.

@param dependNode: The node to get the parent for.
@return: The parent node.
*/
{

	// Verify this is a dag node
	//
	if (!dependNode.hasFn(MFn::kDagNode))
	{

		return MObject::kNullObj;

	}

	// Initialize function set
	//
	MDagPath dagPath = TransformConstraint::getAPathTo(dependNode);
	MFnDagNode fnDagNode(dagPath);

	unsigned int parentCount = fnDagNode.parentCount();

	if (parentCount == 1)
	{

		return fnDagNode.parent(0);

	}
	else
	{

		return MObject::kNullObj;

	}

}


MObject TransformConstraint::getNodeByUUID(MUuid uuid)
{

	// Collect all nodes from the given UUID
	//
	MObjectArray dependNodes = TransformConstraint::getNodesByUUID(uuid);

	// Check if this node has been referenced
	//
	MObject dependNode = this->thisMObject();
	MFnDependencyNode fnDependNode(dependNode);

	if (fnDependNode.isFromReferencedFile())
	{

		// Find reference associated with this node
		//
		MObject reference = TransformConstraint::getAssociatedReferenceNode(dependNode);
		unsigned int numNodes = dependNodes.length();

		switch (numNodes)
		{

			case 0:
				return MObject::kNullObj;
				break;

			case 1:
				return dependNodes[0];
				break;

			default:
				return TransformConstraint::findReferencedNode(reference, dependNodes);
				break;

		}

	}
	else
	{

		// Inspect the number of found nodes
		// This should be unique but accidents can happen...
		//
		unsigned int numNodes = dependNodes.length();

		switch (numNodes)
		{

			case 0:
				return MObject::kNullObj;
				break;

			default:
				return dependNodes[0];
				break;

		}

	}

	return MObject::kNullObj;

};


MObjectArray TransformConstraint::getNodesByUUID(MUuid uuid)
/**
Returns a list of nodes associated with the given UUID

@param uuid: The UUID to test against.
@return: An array of dependency nodes with the given UUID.
*/
{

	MStatus status;

	// Add UUID to selection list
	//
	MSelectionList selection;

	status = selection.add(uuid);
	CHECK_MSTATUS(status);

	// Add nodes to object array
	//
	unsigned int numFound = selection.length();

	MObjectArray dependNodes = MObjectArray(numFound);
	MObject dependNode;

	for (unsigned int i = 0; i < numFound; i++)
	{

		status = selection.getDependNode(i, dependNode);
		CHECK_MSTATUS(status);

		dependNodes[i] = dependNode;

	}

	return dependNodes;

};


MObject TransformConstraint::getAssociatedReferenceNode(MObject& dependNode)
/**
Returns the reference node associated with the given dependency node.
If the dependency node is not referenced then a null object is returned instead!

@param dependNode: The node to find the reference for.
@return: The associated reference node.
*/
{

	MFnDependencyNode fnDependNode(dependNode);

	if (fnDependNode.isFromReferencedFile())
	{

		MItDependencyNodes iterNodes = MItDependencyNodes(MFn::kReference);

		MObject reference;
		MFnReference fnReference;

		while (!iterNodes.isDone())
		{

			reference = iterNodes.thisNode();
			fnReference.setObject(reference);

			if (fnReference.containsNodeExactly(dependNode))
			{

				return reference;

			}

			iterNodes.next();

		}

	}

	return MObject::kNullObj;

};


MObject TransformConstraint::findReferencedNode(MObject& reference, MObjectArray& dependNodes)
/**
Returns the node from the supplied array that belongs to the given reference.

@param reference: The reference node to compare against.
@param dependNodes: An array of dependency nodes to test against.
@return: The node associated with the given reference.
*/
{

	// Iterate through dependency nodes
	//
	MFnReference fnReference(reference);
	unsigned int numDependNodes = dependNodes.length();

	MObject dependNode;

	for (unsigned int i = 0; i < numDependNodes; i++)
	{

		// Check if node belongs to reference
		//
		dependNode = dependNodes[i];

		if (fnReference.containsNodeExactly(dependNode))
		{

			return dependNode;

		}

	}

	return MObject::kNullObj;

}


MMatrix TransformConstraint::createPositionMatrix(MVector position)
/**
Creates a position matrix from the given vector.

@param position: The vector to convert.
@return: The new position matrix.
*/
{

	double matrix[4][4] = {
		{ 1.0, 0.0, 0.0, 0.0 },
		{ 0.0, 1.0, 0.0, 0.0 },
		{ 0.0, 0.0, 1.0, 0.0 },
		{ position.x, position.y, position.z, 1.0 },
	};

	return MMatrix(matrix);

};


MMatrix TransformConstraint::createRotationMatrix(MVector rotation, int rotateOrder)
/**
Creates a rotation matrix from the given vector and rotation order.

@param rotation: The vector to convert.
@param rotateOrder: The order of rotations.
@return: The new rotation matrix.
*/
{

	return MEulerRotation(rotation, MEulerRotation::RotationOrder(rotateOrder)).asMatrix();

};


MMatrix TransformConstraint::createScaleMatrix(double3 scale)
/**
Creates a scale matrix from the given vector.

@param scale: The vector to convert.
@return: The new scale matrix.
*/
{

	double matrix[4][4] = {
		{ scale[0], 0.0, 0.0, 0.0 },
		{ 0.0, scale[1], 0.0, 0.0 },
		{ 0.0, 0.0, scale[2], 0.0 },
		{ 0.0, 0.0, 0.0, 1.0 },
	};

	return MMatrix(matrix);

};


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


MMatrix TransformConstraint::blendMatrices(MMatrix startMatrix, MMatrix endMatrix, float weight)
/**
Interpolates the two given matrices using the supplied weight.
Both translate and scale will be lerp'd while rotation will be slerp'd.

@param startMatrix: The start matrix.
@param endMatrix: The end matrix.
@param weight: The amount to blend.
@return: The interpolated matrix.
*/
{

	// Initialize transformation matrices
	//
	MTransformationMatrix startTransformationMatrix = MTransformationMatrix(startMatrix);
	MTransformationMatrix endTransformationMatrix = MTransformationMatrix(endMatrix);

	// Interpolate translation
	//
	MVector startTranslation = startTransformationMatrix.getTranslation(MSpace::kTransform);
	MVector endTranslation = endTransformationMatrix.getTranslation(MSpace::kTransform);

	MVector translation = lerp(startTranslation, endTranslation, weight);

	// Interpolate rotation
	//
	MQuaternion startQuat = startTransformationMatrix.rotation();
	MQuaternion endQuat = endTransformationMatrix.rotation();

	MQuaternion quat = TransformConstraint::slerp(startQuat, endQuat, weight);

	// Interpolate scale
	//
	double startScale[3];
	double endScale[3];

	startTransformationMatrix.getScale(startScale, MSpace::kTransform);
	endTransformationMatrix.getScale(endScale, MSpace::kTransform);

	double scale[3] = {
		lerp(startScale[0], endScale[0], weight),
		lerp(startScale[1], endScale[1], weight),
		lerp(startScale[2], endScale[2], weight)
	};

	// Compose interpolated matrix
	//
	MMatrix translateMatrix = TransformConstraint::createPositionMatrix(translation);
	MMatrix rotateMatrix = quat.asMatrix();
	MMatrix scaleMatrix = TransformConstraint::createScaleMatrix(scale);

	return scaleMatrix * rotateMatrix * translateMatrix;

};


MMatrix TransformConstraint::blendMatrices(MMatrix restMatrix, MMatrixArray matrices, MFloatArray weights)
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
			MFloatArray clampedWeights = TransformConstraint::clamp(weights);
			float weightSum = TransformConstraint::sum(clampedWeights);

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

				matrix = TransformConstraint::blendMatrices(matrix, matrices[i], weights[i]);

			}

			return matrix;

		}
		break;

	}

	return MMatrix::identity;

};


MQuaternion TransformConstraint::slerp(MQuaternion startQuat, MQuaternion endQuat, float weight)
/**
Spherical interpolates two quaternions.

@param startQuat: Start Quaternion.
@param endQuat: End Quaternion.
@param weight: The amount to interpolate.
@return: The interpolated quaternion.
*/
{

	// Calculate angle between quats
	// If startQuat == endQuat or startQuat == -endQuat then theta = 0 and we can return startQuat
	//
	double cos_half_theta = startQuat.w * endQuat.w + startQuat.x * endQuat.x + startQuat.y * endQuat.y + startQuat.z * endQuat.z;

	if (abs(cos_half_theta) >= 1.0)
	{

		return MQuaternion(startQuat);

	}

	// Calculate temporary values
	// If theta = 180 degrees then result is not fully defined
	// We could rotate around any axis normal to startQuat or endQuat
	//
	MQuaternion quat = MQuaternion();

	double half_theta = acos(cos_half_theta);
	double sin_half_theta = sqrt(1.0 - cos_half_theta * cos_half_theta);

	if (fabs(sin_half_theta) < 0.001)
	{

		quat.x = startQuat.x * 0.5 + endQuat.x * 0.5;
		quat.y = startQuat.y * 0.5 + endQuat.y * 0.5;
		quat.z = startQuat.z * 0.5 + endQuat.z * 0.5;
		quat.w = startQuat.w * 0.5 + endQuat.w * 0.5;

		return quat;

	}

	// Calculate quaternion
	//
	double ratio_a = sin((1.0 - weight) * half_theta) / sin_half_theta;
	double ratio_b = sin(weight * half_theta) / sin_half_theta;

	quat.w = startQuat.w * ratio_a + endQuat.w * ratio_b;
	quat.x = startQuat.x * ratio_a + endQuat.x * ratio_b;
	quat.y = startQuat.y * ratio_a + endQuat.y * ratio_b;
	quat.z = startQuat.z * ratio_a + endQuat.z * ratio_b;

	return quat;

}


float TransformConstraint::sum(MFloatArray items)
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


MFloatArray TransformConstraint::clamp(MFloatArray items)
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


void TransformConstraint::postConstructor()
/**
Internally maya creates two objects when a user defined node is created, the internal MObject and the user derived object.
The association between the these two objects is not made until after the MPxNode constructor is called.
This implies that no MPxNode member function can be called from the MPxNode constructor.
The postConstructor will get called immediately after the constructor when it is safe to call any MPxNode member function.

@return: void
*/
{

	// Store reference to this instance
	//
	this->instances.insert(std::make_pair(this->hashCode(), this));

};


long TransformConstraint::hashCode()
/**
Returns the hash code for this instance.

@return: Hash code.
*/
{

	return MObjectHandle(this->thisMObject()).hashCode();

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
	MFnMessageAttribute fnMessageAttr;

	// Input attributes:
	// ".restTranslateX" attribute
	//
	TransformConstraint::restTranslateX = fnUnitAttr.create("restTranslateX", "rtx", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	// ".restTranslateY" attribute
	//
	TransformConstraint::restTranslateY = fnUnitAttr.create("restTranslateY", "rty", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	// ".restTranslateZ" attribute
	//
	TransformConstraint::restTranslateZ = fnUnitAttr.create("restTranslateZ", "rtz", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	// ".restTranslate" attribute
	//
	TransformConstraint::restTranslate = fnNumericAttr.create("restTranslate", "rt", TransformConstraint::restTranslateX, TransformConstraint::restTranslateY, TransformConstraint::restTranslateZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	// ".restRotateX" attribute
	//
	TransformConstraint::restRotateX = fnUnitAttr.create("restRotateX", "rrx", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	// ".restRotateY" attribute
	//
	TransformConstraint::restRotateY = fnUnitAttr.create("restRotateY", "rry", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	// ".restRotateZ" attribute
	//
	TransformConstraint::restRotateZ = fnUnitAttr.create("restRotateZ", "rrz", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	// ".restRotate" attribute
	//
	TransformConstraint::restRotate = fnNumericAttr.create("restRotate", "rr", TransformConstraint::restRotateX, TransformConstraint::restRotateY, TransformConstraint::restRotateZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	// ".restScaleX" attribute
	//
	TransformConstraint::restScaleX = fnNumericAttr.create("restScaleX", "rsx", MFnNumericData::kDouble, 1.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	// ".restScaleY" attribute
	//
	TransformConstraint::restScaleY = fnNumericAttr.create("restScaleY", "rsy", MFnNumericData::kDouble, 1.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	// ".restScaleZ" attribute
	//
	TransformConstraint::restScaleZ = fnNumericAttr.create("restScaleZ", "rsz", MFnNumericData::kDouble, 1.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	// ".restScale" attribute
	//
	TransformConstraint::restScale = fnNumericAttr.create("restScale", "rs", TransformConstraint::restScaleX, TransformConstraint::restScaleY, TransformConstraint::restScaleZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	// ".constraintObject" attribute
	//
	TransformConstraint::constraintObject = fnMessageAttr.create("constraintObject", "co", &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

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

	// ".constraintJointOrientX" attribute
	//
	TransformConstraint::constraintJointOrientX = fnUnitAttr.create("constraintJointOrientX", "cjox", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	// ".constraintJointOrientY" attribute
	//
	TransformConstraint::constraintJointOrientY = fnUnitAttr.create("constraintJointOrientY", "cjoy", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	// ".constraintJointOrientZ" attribute
	//
	TransformConstraint::constraintJointOrientZ = fnUnitAttr.create("constraintJointOrientZ", "cjoz", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	// ".constraintJointOrient" attribute
	//
	TransformConstraint::constraintJointOrient = fnNumericAttr.create("constraintJointOrient", "cjo", TransformConstraint::constraintJointOrientX, TransformConstraint::constraintJointOrientY, TransformConstraint::constraintJointOrientZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	// ".constraintParentInverseMatrix" attribute
	//
	TransformConstraint::constraintParentInverseMatrix = fnMatrixAttr.create("constraintParentInverseMatrix", "cpim", MFnMatrixAttribute::kDouble, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	// Target Attributes:
	// Define ".targetWeight" attribute
	//
	TransformConstraint::targetWeight = fnNumericAttr.create("targetWeight", "tw", MFnNumericData::kFloat, 0.0f, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.setMin(0.0));
	CHECK_MSTATUS(fnNumericAttr.setMax(1.0));

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetTranslateX" attribute
	//
	TransformConstraint::targetTranslateX = fnUnitAttr.create("targetTranslateX", "ttx", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetTranslateY" attribute
	//
	TransformConstraint::targetTranslateY = fnUnitAttr.create("targetTranslateY", "tty", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetTranslateZ" attribute
	//
	TransformConstraint::targetTranslateZ = fnUnitAttr.create("targetTranslateZ", "ttz", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetTranslate" attribute
	//
	TransformConstraint::targetTranslate = fnNumericAttr.create("targetTranslate", "tt", TransformConstraint::targetTranslateX, TransformConstraint::targetTranslateY, TransformConstraint::targetTranslateZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetOffsetTranslateX" attribute
	//
	TransformConstraint::targetOffsetTranslateX = fnUnitAttr.create("targetOffsetTranslateX", "totx", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetOffsetTranslateY" attribute
	//
	TransformConstraint::targetOffsetTranslateY = fnUnitAttr.create("targetOffsetTranslateY", "toty", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetOffsetTranslateZ" attribute
	//
	TransformConstraint::targetOffsetTranslateZ = fnUnitAttr.create("targetOffsetTranslateZ", "totz", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetOffsetTranslate" attribute
	//
	TransformConstraint::targetOffsetTranslate = fnNumericAttr.create("targetOffsetTranslate", "tot", TransformConstraint::targetOffsetTranslateX, TransformConstraint::targetOffsetTranslateY, TransformConstraint::targetOffsetTranslateZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetJointOrientX" attribute
	//
	TransformConstraint::targetJointOrientX = fnUnitAttr.create("targetJointOrientX", "tjox", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetJointOrientY" attribute
	//
	TransformConstraint::targetJointOrientY = fnUnitAttr.create("targetJointOrientY", "tjoy", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetJointOrientZ" attribute
	//
	TransformConstraint::targetJointOrientZ = fnUnitAttr.create("targetJointOrientZ", "tjoz", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetJointOrient" attribute
	//
	TransformConstraint::targetJointOrient = fnNumericAttr.create("targetJointOrient", "tjo", TransformConstraint::targetJointOrientX, TransformConstraint::targetJointOrientY, TransformConstraint::targetJointOrientZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetRotateX" attribute
	//
	TransformConstraint::targetRotateX = fnUnitAttr.create("targetRotateX", "trx", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetRotateY" attribute
	//
	TransformConstraint::targetRotateY = fnUnitAttr.create("targetRotateY", "try", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetRotateZ" attribute
	//
	TransformConstraint::targetRotateZ = fnUnitAttr.create("targetRotateZ", "trz", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetRotate" attribute
	//
	TransformConstraint::targetRotate = fnNumericAttr.create("targetRotate", "tr", TransformConstraint::targetRotateX, TransformConstraint::targetRotateY, TransformConstraint::targetRotateZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetOffsetRotateX" attribute
	//
	TransformConstraint::targetOffsetRotateX = fnUnitAttr.create("targetOffsetRotateX", "torx", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetOffsetRotateY" attribute
	//
	TransformConstraint::targetOffsetRotateY = fnUnitAttr.create("targetOffsetRotateY", "tory", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetOffsetRotateZ" attribute
	//
	TransformConstraint::targetOffsetRotateZ = fnUnitAttr.create("targetOffsetRotateZ", "torz", MFnUnitAttribute::kAngle, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetOffsetRotate" attribute
	//
	TransformConstraint::targetOffsetRotate = fnNumericAttr.create("targetOffsetRotate", "tor", TransformConstraint::targetOffsetRotateX, TransformConstraint::targetOffsetRotateY, TransformConstraint::targetOffsetRotateZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

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

	CHECK_MSTATUS(fnEnumAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetScaleX" attribute
	//
	TransformConstraint::targetScaleX = fnNumericAttr.create("targetScaleX", "tsx", MFnNumericData::kDouble, 1.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetScaleY" attribute
	//
	TransformConstraint::targetScaleY = fnNumericAttr.create("targetScaleY", "tsy", MFnNumericData::kDouble, 1.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetScaleZ" attribute
	//
	TransformConstraint::targetScaleZ = fnNumericAttr.create("targetScaleZ", "tsz", MFnNumericData::kDouble, 1.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetScale" attribute
	//
	TransformConstraint::targetScale = fnNumericAttr.create("targetScale", "ts", TransformConstraint::targetScaleX, TransformConstraint::targetScaleY, TransformConstraint::targetScaleZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetOffsetScaleX" attribute
	//
	TransformConstraint::targetOffsetScaleX = fnNumericAttr.create("targetOffsetScaleX", "tosx", MFnNumericData::kDouble, 1.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetOffsetScaleY" attribute
	//
	TransformConstraint::targetOffsetScaleY = fnNumericAttr.create("targetOffsetScaleY", "tosy", MFnNumericData::kDouble, 1.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetOffsetScaleZ" attribute
	//
	TransformConstraint::targetOffsetScaleZ = fnNumericAttr.create("targetOffsetScaleZ", "tosz", MFnNumericData::kDouble, 1.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetOffsetScale" attribute
	//
	TransformConstraint::targetOffsetScale = fnNumericAttr.create("targetOffsetScale", "tos", TransformConstraint::targetOffsetScaleX, TransformConstraint::targetOffsetScaleY, TransformConstraint::targetOffsetScaleZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetRotatePivotX" attribute
	//
	TransformConstraint::targetRotatePivotX = fnUnitAttr.create("targetRotatePivotX", "trpx", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetRotatePivotY" attribute
	//
	TransformConstraint::targetRotatePivotY = fnUnitAttr.create("targetRotatePivotY", "trpy", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetRotatePivotZ" attribute
	//
	TransformConstraint::targetRotatePivotZ = fnUnitAttr.create("targetRotatePivotZ", "trpz", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetRotatePivot" attribute
	//
	TransformConstraint::targetRotatePivot = fnNumericAttr.create("targetRotatePivot", "trp", TransformConstraint::targetRotatePivotX, TransformConstraint::targetRotatePivotY, TransformConstraint::targetRotatePivotZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetRotateTranslateX" attribute
	//
	TransformConstraint::targetRotateTranslateX = fnUnitAttr.create("targetRotateTranslateX", "trptx", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetRotateTranslateY" attribute
	//
	TransformConstraint::targetRotateTranslateY = fnUnitAttr.create("targetRotateTranslateY", "trpty", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetRotateTranslateZ" attribute
	//
	TransformConstraint::targetRotateTranslateZ = fnUnitAttr.create("targetRotateTranslateZ", "trptz", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetRotateTranslate" attribute
	//
	TransformConstraint::targetRotateTranslate = fnNumericAttr.create("targetRotateTranslate", "trpt", TransformConstraint::targetRotateTranslateX, TransformConstraint::targetRotateTranslateY, TransformConstraint::targetRotateTranslateZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetScalePivotX" attribute
	//
	TransformConstraint::targetScalePivotX = fnUnitAttr.create("targetScalePivotX", "tspx", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetScalePivotY" attribute
	//
	TransformConstraint::targetScalePivotY = fnUnitAttr.create("targetScalePivotY", "tspy", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetScalePivotZ" attribute
	//
	TransformConstraint::targetScalePivotZ = fnUnitAttr.create("targetScalePivotZ", "tspz", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetScalePivot" attribute
	//
	TransformConstraint::targetScalePivot = fnNumericAttr.create("targetScalePivot", "tsp", TransformConstraint::targetScalePivotX, TransformConstraint::targetScalePivotY, TransformConstraint::targetScalePivotZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetScaleTranslateX" attribute
	//
	TransformConstraint::targetScaleTranslateX = fnUnitAttr.create("targetScaleTranslateX", "tsptx", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetScaleTranslateY" attribute
	//
	TransformConstraint::targetScaleTranslateY = fnUnitAttr.create("targetScaleTranslateY", "tspty", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetScaleTranslateZ" attribute
	//
	TransformConstraint::targetScaleTranslateZ = fnUnitAttr.create("targetScaleTranslateZ", "tsptz", MFnUnitAttribute::kDistance, 0.0, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnUnitAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetScaleTranslate" attribute
	//
	TransformConstraint::targetScaleTranslate = fnNumericAttr.create("targetScaleTranslate", "tspt", TransformConstraint::targetScaleTranslateX, TransformConstraint::targetScaleTranslateY, TransformConstraint::targetScaleTranslateZ, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	CHECK_MSTATUS(fnNumericAttr.addToCategory(TransformConstraint::targetCategory));

	// ".targetParentMatrix" attribute
	//
	TransformConstraint::targetParentMatrix = fnMatrixAttr.create("targetParentMatrix", "tpm", MFnMatrixAttribute::kDouble, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

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
	CHECK_MSTATUS(fnCompoundAttr.addChild(TransformConstraint::targetRotateTranslate));
	CHECK_MSTATUS(fnCompoundAttr.addChild(TransformConstraint::targetScalePivot));
	CHECK_MSTATUS(fnCompoundAttr.addChild(TransformConstraint::targetScaleTranslate));
	CHECK_MSTATUS(fnCompoundAttr.addChild(TransformConstraint::targetParentMatrix));

	CHECK_MSTATUS(fnCompoundAttr.setArray(true));

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
	CHECK_MSTATUS(TransformConstraint::addAttribute(TransformConstraint::constraintObject));

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

	// Create child added callback
	//
	TransformConstraint::childAddedCallbackId = MDagMessage::addChildAddedCallback(onChildAdded);

	return MS::kSuccess;

};
