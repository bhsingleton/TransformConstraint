//
// File: pluginMain.cpp
//
// Author: Ben Singleton
//

#include "TransformConstraint.h"

#include <maya/MFnPlugin.h>


MStatus initializePlugin(MObject obj)
{ 

	MStatus status;

	MFnPlugin plugin(obj, "Ben Singleton", "2020", "Any");
	status = plugin.registerNode("transformConstraint", TransformConstraint::id, TransformConstraint::creator, TransformConstraint::initialize, MPxNode::kConstraintNode);
	
	if (!status)
	{

		status.perror("registerNode");
		return status;

	}

	return status;

}


MStatus uninitializePlugin(MObject obj)
{

	MStatus status;

	MFnPlugin plugin(obj);
	status = plugin.deregisterNode(TransformConstraint::id);

	if (!status)
	{

		status.perror("deregisterNode");
		return status;

	}

	status = MDagMessage::removeCallback(TransformConstraint::childAddedCallbackId);

	if (!status)
	{

		status.perror("removeCallback");
		return status;

	}

	return status;

}
