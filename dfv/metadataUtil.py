#!/bin/env python
# coding: utf-8

from collections import OrderedDict, namedtuple
import re
import copy
import os.path

class MetadataField():

    def __init__(self, name, description, qualifier, feature, entry, type_="string", mss_required=False, dfast_required=False, pattern=re.compile(r".*"), value=""):
        """
            name: a name for a variable used in DFAST (case sensitive, must be unique)
            description: description of this field
                ex) name: locusTagPrefix  -> description: Locus Tag Prefix
            qualifier is a term used in MSS format. ex) ab_name
            feature: a feature label where this metadata are described. ex) SUBMITTER, REFERENCE, ST_COMMENT
            entry: an entry label where this metadata are described.  COMMON|SEQUENCE|OTHER
            array: Is set True if multiple values are acceptable. Must be splitted with semicolons. ex) "spam; ham; eggs"
            mss_required: Is used to check if a generated file is in the valid DDBJ-MSS format.
            dfast_required: Is set True if a null value is acceptable. Only used for GUI.
            pattern: a regex pattern that must MATCH the value
            value: a value for this field
        """
        if type_ == "array":
            value = [x.strip() for x in value.split(";") if x.strip()]
        self.name, self.description, self.qualifier, self.feature, self.entry, self.type, \
            self.mss_required, self.dfast_required, self.pattern, self.value = \
            name, description, qualifier, feature, entry, type_, mss_required, dfast_required, pattern, value

    def __repr__(self):
        if self.type == "array":
            values = "; ".join(self.value)
            return "{self.name}\t{values}".format(self=self, values=values)
        else:
            return "{self.name}\t{self.value}".format(self=self)

    def appendValue(self, value):
        assert self.type == "array"
        self.value.append(value)

    def setValue(self, value):
        if self.type == "array":
            self.value = [x.strip() for x in value.split(";") if x.strip()]
        else:
            self.value = value

    def getAdditionalField(self, newName, newValue=""):
        if ":" not in newName:
            print("no colon : in the parameter name")
            raise AssertionError
        num = newName.split(":")[-1]
        newFeature = "{0}:{1}".format(self.feature, newName.split(":")[-1])
        newObj = MetadataField(newName, self.description, self.qualifier, newFeature, self.entry,
                               self.type, self.mss_required, self.dfast_required, self.pattern, newValue)
        return newObj

    def validate(self, ignoreNull=True):
        # if ignoreNull and (not self.mss_required) and (self.value == "" or self.value == []):
        if self.value == "" or self.value == []:
            if ignoreNull or (not self.mss_required):
                return {}
            else:
                # return {self.name: {"msg": "Empty value", "qualifier": self.qualifier, "feature": self.feature}}
                return {"MISSING_VALUES": {"key": self.name}}

        failedValues = {}
        if self.type == "array":
            for eachVal in self.value:
                if self.pattern.match(eachVal) is None:
                    failedValues.setdefault("values", []).append(eachVal)
                    failedValues["pattern"] = self.pattern.pattern.replace("^(", "").replace(")$", "")

        else:
            if self.pattern.match(self.value) is None:
                failedValues.setdefault("values", []).append(self.value)
                failedValues["pattern"] = self.pattern.pattern.replace("^(", "").replace(")$", "")

        if len(failedValues) > 0:
            failedValues.update({"key": self.name})
            return {"INCORRECT_VALUES": failedValues}
            # failedValues.update({"qualifier": self.qualifier, "feature": self.feature})
            # return {self.name: failedValues}
        else:
            return {}

    def render(self, addFeature=False):
        # TODO: cuttently addFeature is always False. It can be removed in the future version.
        RET = []
        retFeature = self.feature if addFeature else ""

        if self.value == "" or self.value == []:
            RET.append(["", retFeature, "", self.qualifier, ""])
            if self.qualifier == "ab_name":
                # if ab_name has no value, render additional 2 empty lines into the template
                RET.append(["", "", "", self.qualifier, ""])
                RET.append(["", "", "", self.qualifier, ""])
        elif self.type == "array":
            if self.feature == "ST_COMMENT":
                RET.append(["", "", "", self.qualifier, "; ".join(self.value)])
            else:
                for eachVal in self.value:
                    if eachVal:
                        RET.append(["", "", "", self.qualifier, eachVal])
            RET[0][1] = retFeature
        elif self.type == "boolean" and self.value:
            if self.value != "NO":
                RET.append(["", retFeature, "", self.qualifier, ""])            
        else:
            RET.append(["", retFeature, "", self.qualifier, self.value])
        return RET

class Metadata():
    METADATA_DEFINITION_FILE = os.path.dirname(os.path.abspath(__file__)) + "/metadataDefinition.tsv"
    METADATA_RULES_FILE = os.path.dirname(os.path.abspath(__file__)) + "/metadataRules.tsv"

    def __init__(self):
        self.fields = OrderedDict()
        self.refNumber = 0
        self.commentNumber = 0

    @staticmethod
    def initialize():
        metadata = Metadata()
        for line in open(metadata.__class__.METADATA_DEFINITION_FILE):
            if line.startswith("#"):
                continue
            else:
                name, description, qualifier, feature, entry, type_, mss_required, dfast_required, pattern, value = line.strip("\n").split("\t")[0:10]
                if entry == "EX_SOURCE":
                    continue  # skipe extended qualifiers for source features
                # if array == "array":
                #     array = True
                #     # value = [x.strip() for x in value.split(";") if x.strip()]
                # else:
                #     array = False
                mss_required = True if mss_required == "TRUE" else False
                dfast_required = True if dfast_required == "TRUE" else False
                pattern = re.compile(r"^({0})$".format(pattern))
                field = MetadataField(name, description, qualifier, feature, entry, type_, mss_required, dfast_required, pattern, value)
                metadata.fields[name] = field
        return metadata

    @staticmethod
    def getDefaultDefinition(projectType="other"):
        '''for mtg only'''
        metadata = Metadata.initializeFromFormData({"projectType": projectType})
        return metadata.getDefinition()


    def getDefinition(self):
        ret = []
        for field in self.fields.values():
            if field.entry == "DFAST" and field.name != "projectType":

                continue
            pattern = field.pattern.pattern.replace("^(", "").replace(")$", "")

            ret.append({"name": field.name, "description": field.description, "qualifier": field.qualifier, "feature": field.feature,
                         "entry": field.entry, "type": field.type, "mss_required": field.mss_required, 
                         "dfast_required": field.dfast_required, "pattern": pattern, "value":field.value })
        sourceQualifiers = self.getSourceQualifiers()
        return ret, sourceQualifiers

    def getSourceQualifiers(self):
        ret = []
        for line in open(Metadata.METADATA_DEFINITION_FILE):
            if line.startswith("#"):
                continue
            else:
                name, description, qualifier, feature, entry, type_, mss_required, dfast_required, pattern, value = line.strip("\n").split("\t")[0:10]
                if feature == "source" and name not in self.fields:
                    ret.append({"name":name, "description":description})
        return ret

    # @staticmethod
    # def sourceQualifiers():
    #     ret = []
    #     # metadata = Metadata()
    #     for line in open(Metadata.METADATA_DEFINITION_FILE):
    #         if line.startswith("#"):
    #             continue
    #         else:
    #             name, description, qualifier, feature, entry, type_, mss_required, dfast_required, pattern, value = line.strip("\n").split("\t")[0:10]
    #             if feature != "source":
    #                 continue  # skipe extended qualifiers for source features

    #             # mss_required = True if mss_required == "TRUE" else False
    #             # dfast_required = True if dfast_required == "TRUE" else False
    #             # pattern = re.compile(r"^({0})$".format(pattern))
    #             # field = MetadataField(name, description, qualifier, feature, entry, type_, mss_required, dfast_required, pattern, value)
    #             # metadata.fields[name] = field
    #             ret.append(name)
    #     return ret

    # @staticmethod
    # def definition_old():
    #     ret = []
    #     for line in open(Metadata.METADATA_DEFINITION_FILE):
    #         if line.startswith("#"):
    #             continue
    #         else:
    #             name, description, qualifier, feature, entry, type_, mss_required, dfast_required, pattern, value = line.strip("\n").split("\t")[0:10]
    #             if entry == "EX_SOURCE":
    #                 continue  # skipe extended qualifiers for source features
    #             mss_required = True if mss_required == "TRUE" else False
    #             dfast_required = True if dfast_required == "TRUE" else False

    #             ret.append({"name": name, "description": description, "qualifier": qualifier,
    #                          "entry": entry, "type": type_, "mss_required": mss_required, 
    #                          "dfast_required": dfast_required, "pattern": pattern})
    #     return ret

    def getValue(self, fieldName, defaultVal="", idx=0):
        field = self.fields.get(fieldName)
        if field is None:
            return defaultVal
        if field.type == "array":
            if len(field.value) > idx:
                return field.value[idx]
            else:
                return defaultVal
        else:
            return field.value or defaultVal

    def getValues(self, fieldName):
        field = self.fields.get(fieldName)
        if field is None or field.type != "array":
            return []
        else:
            return field.value

    def setValue(self, fieldName, val):
        field = self.fields.get(fieldName)
        if field:
            field.setValue(val)
        else:
            raise IndexError("Field name not found")

    def setDefaultValue(self, fieldName, val):
        '''if value is empty, default value will be set'''
        if self.getValue(fieldName) == "":
            self.setValue(fieldName, val)

    def addField(self, name, value=None, mss_required=None, dfast_required=None, pattern=None):
        '''Add field and set value if field does not exist.
        If field already exists, only its value is set'''
        if name in self.fields:
            if value is not None:
                self.setDefaultValue(name, value)
        else:
            for line in open(self.__class__.METADATA_DEFINITION_FILE):
                cols = line.strip("\n").split("\t")
                if name == cols[0]:
                    name, description, qualifier, feature, entry, type_, default_mss_required, default_dfast_required, default_pattern, default_value = cols[0:10]

                    mss_required = default_mss_required if mss_required is None else mss_required
                    dfast_required = default_dfast_required if dfast_required is None else dfast_required
                    pattern = default_pattern if pattern is None else pattern
                    pattern = re.compile(r"^({0})$".format(pattern))

                    mss_required = True if mss_required == "TRUE" else False
                    dfast_required = True if dfast_required == "TRUE" else False

                    field = MetadataField(name, description, qualifier, feature, entry, type_, mss_required, dfast_required, pattern, default_value)
                    if value is not None:
                        field.setValue(value)

                    self.fields[name] = field


    def getFields(self):
        return list(self.fields.keys())  # self.fields is a dictionary of metadata

    def toTSV(self, outputFileName):
        Buffer = ""
        for field in self.fields.values():
            Buffer += repr(field) + "\n"
        with open(outputFileName, "w") as f:
            f.write(Buffer)

    def validateValues(self, ignoreNull=True):

        def checkSeqNumberInconsistency():
            if self.getValue("projectType") == "gnm":
                if not (len(self.fields["seqNames"].value) == len(self.fields["seqTopologies"].value) == len(self.fields["seqTypes"].value)):
                    return {"message": "seqNames, seqTopologies, seqTypes must be the same length.", "keys": ["seqNames", "seqTopologies", "seqTypes"]}

        errors = {"status": "success",
                  "MISSING_VALUES": [],
                  "INCORRECT_VALUES": [],
                  "INCONSISTENT_VALUES": [],
                  }
        for field in self.fields.values():
            ret = field.validate(ignoreNull)
            # ret should be like: {"INCORRECT_VALUE": {"key": hoge, "value": [ham, spam, egg]}
            for key, val in ret.items():
                errors[key].append(val)

        err = checkSeqNumberInconsistency()
        if err:
            errors["INCONSISTENT_VALUES"].append(err)

        if len(errors["MISSING_VALUES"]) > 0 or len(errors["INCORRECT_VALUES"]) > 0 or len(errors["INCONSISTENT_VALUES"]) > 0:
            errors["status"] = "fail"
        return errors
        # failedFields = OrderedDict()
        # for field in self.fields.values():
        #     ret = field.validate(ignoreNull)
        #     for fieldName, result in ret.items():
        #         failedFields[fieldName] = result
        # for key, value in failedFields.items():
        # print(key, value)
        # return failedFields

    def renderSourceFeature(self, seqNum=-1, seqRecord=None):
        def parseSeqRecord(seqRecord):
            if seqRecord is None or seqNum < 0:
                return "", "E"
            else:
                seqName = seqRecord.name
                length = len(seqRecord) if seqRecord else "E"
                return seqName, length

        def renderField_if_hasValue(outputData, fieldName):
            field = self.fields.get(fieldName)
            if field and not (field.value == "" or field.value == []):
                outputData += field.render()
            elif field and field.type == "boolean" and field.value:
                outputData += field.render()

            elif field and field.mss_required:
                outputData += field.render()

        entryName, length = parseSeqRecord(seqRecord)
        seqStatus = self.getValue("seqStatus", "draft")
        organismName = self.getValue("organism", "")
        if not organismName and (self.getValue("genus", "") and self.getValue("species", "")):
            organismName = self.getValue("genus", "Genus") + " " + self.getValue("species", "sp.")
        # organismName = self.getValue("organism", "") or self.getValue("genus", "Genus") + " " + self.getValue("species", "sp.")

        if seqStatus == "complete":
            if seqNum >= 0:
                seqTopology = self.getValue("seqTopologies", "linear", idx=seqNum)
                seqType = self.getValue("seqTypes", "unkown", idx=seqNum)
                # seqName = self.getValue("seqNames", "unkown", idx=seqNum)
            else:
                seqTopology, seqType = "linear", ""
        else:
            seqTopology, seqType = "linear", ""
            seqRank = self.getValue("seqRank", "contig")

        outputData = []
        outputData.append([entryName, "source", "1..{}".format(length), "organism", organismName])
        renderField_if_hasValue(outputData, "strain")
        renderField_if_hasValue(outputData, "isolate")
        renderField_if_hasValue(outputData, "cultivar")
        renderField_if_hasValue(outputData, "mol_type")

        if self.getValue("typeStrain", "NO") == "YES":
            # outputData.append(["", "", "", "type_material", "type strain of " + organismName])
            outputData.append(["", "", "", "type_material", "type strain"])

        # renderField_if_hasValue(outputData, "cultureCollection")

        if seqType == "plasmid":
            outputData.append(["", "", "", "plasmid", entryName])
        elif seqType == "unplaced":
            outputData.append(["", "", "", "note", "unplaced contig"])

        ff_definition = self.getValue("ff_definition", "")
        projectType = self.getValue("projectType", "wgs")
        if projectType == "gnm":
            pass
            # if ff_definition == "":
            #     # ff_definition = "@@[organism]@@ DNA, strain: @@[strain]@@, complete genome, @@[entry]@@"
            #     ff_definition = "@@[organism]@@ @@[strain]@@ DNA, complete genome: @@[entry]@@"
            # outputData.append(["", "", "", "ff_definition", ff_definition])
        elif projectType in ["wgs", "tsa"] :  # in case of "draft":
            # ff_definition = "@@[organism]@@ @@[strain]@@ DNA, @@[submitter_seqid]@@"
            # if ff_definition == "":
            #     ff_definition = "@@[organism]@@ DNA, strain: @@[strain]@@, %s: @@[entry]@@" % seqRank
            # elif seqRank == "scaffold":
            #     ff_definition = ff_definition.replace("contig", "scaffold")
            # outputData.append(["", "", "", "ff_definition", ff_definition])
            outputData.append(["", "", "", "submitter_seqid", "@@[entry]@@"])
        else:
            pass
            # outputData.append(["", "", "", "ff_definition", ff_definition])
        outputData.append(["", "", "", "ff_definition", ff_definition])

        for field in self.fields.values():
            if field.entry == "EX_SOURCE" and field.name != "typeStrain":
                # print("ADDING EX_SOURCE in renderSource", field.name, field.mss_required, field.value)
                renderField_if_hasValue(outputData, field.name)

        if seqTopology == "circular":
            outputData.append(["", "TOPOLOGY", "", "circular", ""])

        return outputData

    def renderCommonEntry(self, includeSource=False):

        def renderFeature(outputData, featureType):
            # assert featureType in ["DBLINK", "SUBMITTER", "REFERENCE"]
            RET = []
            for field in self.fields.values():
                if field.feature == featureType:
                    if field.value or field.mss_required:
                        RET += field.render()
            if featureType == "ST_COMMENT" and len(RET) == 1:
                # if ST_COMMENT empty, leave without doing anything
                return
            if len(RET) > 0:
                RET[0][1] = featureType.split(":")[0]
                outputData += RET

        outputData = []

        renderFeature(outputData, "DIVISION")
        renderFeature(outputData, "DATATYPE")
        renderFeature(outputData, "KEYWORD")
        renderFeature(outputData, "DBLINK")
        renderFeature(outputData, "SUBMITTER")
        renderFeature(outputData, "REFERENCE")   # render REFERENCE regardless of refNumber
        for i in range(1, self.refNumber):   # render additional REFERENCES
            renderFeature(outputData, "REFERENCE:" + str(i + 1))
        for i in range(self.commentNumber):
            if i == 0:
                renderFeature(outputData, "COMMENT")
            else:
                renderFeature(outputData, "COMMENT:" + str(i + 1))
        renderFeature(outputData, "COMMENT:EST")
        # estcomment = self.fields.get("estcomment")
        # if estcomment and estcomment.value:
        #     outputData += estcomment.render(addFeature=True)

        # print self.getValue("estcomment")
        renderFeature(outputData, "ST_COMMENT")
        renderFeature(outputData, "DATE")

        outputData[0][0] = "COMMON"

        if includeSource:
            outputData += self.renderSourceFeature()
        return outputData
        # debug
        # for line in outputData:
        # print(line)

    def render4DFAST(self, includeSource=False, dfc_annotation=True):
        def appendNewComment(comment):
            if self.commentNumber == 0:
                self.fields["comment"].setValue(comment)
            else:
                fieldName = "comment:{}".format(str(self.commentNumber + 1))
                originalField = self.fields["comment"]
                newCommentField = originalField.getAdditionalField(fieldName, comment)
                self.fields[fieldName] = newCommentField
            self.commentNumber += 1

        def dropComment():
            if self.commentNumber == 0:
                pass
                raise AssertionError("No comment exists.")
            elif self.commentNumber == 1:
                self.fields["comment"].setValue("")
            else:
                del self.fields["comment:{}".format(str(self.commentNumber))]
            self.commentNumber -= 1

        # seqStatus = self.getValue("seqStatus", "draft")
        # if seqStatus == "complete":
        #     self.setValue("division", "")
        #     self.setValue("dataType", "")
        #     self.setValue("keyword", "")
        if dfc_annotation:
            print("xxxxxxx")
            appendNewComment("Annotated by DFAST https://dfast.ddbj.nig.ac.jp/")

        output = self.renderCommonEntry(includeSource)
        dropComment()

        return output

    def addComments(self, metadataDict):
        # print(metadataDict)
        estcomment = metadataDict.get("estcomment") 
        if estcomment:
            self.fields["estcomment"].setValue(estcomment)
        checkNext = True
        while checkNext:
            if self.commentNumber == 0:
                fieldName = "comment"
                value = metadataDict.get("comment") or metadataDict.get("comment:1")

                if value:
                    self.fields[fieldName].setValue(value)
                    self.commentNumber += 1
                    # print("comment 1. added", value)
                else:
                    self.fields[fieldName].value = []
                    self.commentNumber += 1
                    # print("comment 1. not found. exit")

                    # checkNext = False
            else:
                fieldName = "comment:{}".format(str(self.commentNumber + 1))
                value = metadataDict.get(fieldName)
                originalField = self.fields["comment"]
                if fieldName in metadataDict:
                    # print("additional comment. added", value)
                    commentField = originalField.getAdditionalField(fieldName, value)
                    self.fields[fieldName] = commentField
                    self.commentNumber += 1
                else:
                    # print("comment . not found. exit", fieldName, value)
                    checkNext = False

    def addReferences(self, metadataDict):

        def searchPMID(pmid):
            D = OrderedDict()
            D["pubmed"] = pmid
            D["reference"] = "Reference title. PMID:{}".format(pmid)
            D["author"] = "Auhor,{0}-1.; Author,{0}-2.; Author,{0}-3.".format(pmid)
            D["refconsrtm"] = "Consortium dummy"
            D["status"] = "Published"
            D["year"] = "2020"
            D["journal"] = "Journal title. PMID:{}".format(pmid)
            D["volume"] = "V20"
            D["start_page"] = "P100"
            D["end_page"] = "P200"
            return D

        fields = [x for x in self.fields.values() if x.feature == "REFERENCE"]
        checkNext = True
        while checkNext:

            if self.refNumber == 0:  # At least one reference is required, so checkNext will be "True" in this case
                pmid = metadataDict.get("pubmed") or metadataDict.get("pubmed:1")
                if pmid:
                    refInfo = searchPMID(pmid)
                    if "reference" in refInfo:
                        self.refNumber += 1
                        for name, value in refInfo.items():
                            self.fields[name].setValue(value)
                else:
                    for field in fields:
                        value = metadataDict.get(field.name) or metadataDict.get(field.name + ":1", "")
                        self.fields[field.name].setValue(value)
                    self.refNumber += 1
            else:  # If no additional reference is found, so checkNext will be "False"
                pmid = metadataDict.get("pubmed:" + str(self.refNumber + 1))
                if pmid:
                    refInfo = searchPMID(pmid)
                    if "reference" in refInfo:
                        self.refNumber += 1
                        for name, value in refInfo.items():
                            newName = "{}:{}".format(name, str(self.refNumber))
                            additionalField = self.fields[name].getAdditionalField(newName, value)
                            self.fields[newName] = additionalField
                else:
                    if "reference:" + str(self.refNumber + 1) in metadataDict:
                        self.refNumber += 1
                        for field in fields:
                            fieldName = "{}:{}".format(field.name, str(self.refNumber))
                            value = metadataDict.get(fieldName, "")
                            additionalField = field.getAdditionalField(fieldName, value)
                            self.fields[fieldName] = additionalField

                    else:
                        checkNext = False

    # def addReferencesByPMID(self, metadataDict):
    #     def searchPMID(pmid):

    #         return {"reference": "Reference title. PMID:{}".format(pmid),
    #                 "author": "Auhor, {0}-1.; Author, {0}-2.; Author, {0}-3.".format(pmid),
    #                 "refconsrtm": "Consortium dummy",
    #                 "status": "Published",
    #                 "year": "2020",
    #                 "journal": "Journal title. PMID:{}".format(pmid),
    #                 "volume": "V20", "start_page": "P100", "end_page": "P200"}

    #     fields = [x for x in self.fields.values() if x.feature == "REFERENCE"]
    #     checkNext = True
    #     while checkNext:

    #         if self.refNumber == 0:  # At least one reference is required, so checkNext will be "True" in this case
    #             pmid = metadataDict.get("pubmed") or metadataDict.get("pubmed:1")
    #             if pmid:
    #                 refInfo = searchPMID(pmid)
    #                 if "reference" in refInfo:
    #                     self.refNumber += 1
    #                     for name, value in refInfo.items():
    #                         self.fields[name].setValue(value)
    #         else:  # If no additional reference is found, so checkNext will be "False"
    #             pmid = metadataDict.get("pubmed:" + str(self.refNumber + 1))
    #             if pmid:
    #                 refInfo = searchPMID(pmid)
    #                 if "reference" in refInfo:
    #                     self.refNumber += 1
    #                     for name, value in refInfo.values():
    #                         newName = name + str(self.refNumber)
    #                         additionalField = self.fields[name].getAdditionalField(newName, value)
    #                         self.fields[newName] = additionalField
    #                     print("reference added from PUBMED")
    #                 else:
    #                     checkNext = False

    def addExtendedSourceQualifiers(self, formData):
        fields = {}
        # print("DEBUG1:", formData.get("clone"))
        for line in open(self.__class__.METADATA_DEFINITION_FILE):
            if line.startswith("#"):
                continue
            else:
                name, description, qualifier, feature, entry, type_, mss_required, dfast_required, pattern, value = line.strip("\n").split("\t")[0:10]
                if entry != "EX_SOURCE":
                    continue  # skipe non-extended qualifiers for source features
                # if array == "TRUE":
                #     array = True
                    # value = [x.strip() for x in value.split(";") if x.strip()]
                # else:
                #     array = False
                mss_required = True if mss_required == "TRUE" else False
                dfast_required = True if dfast_required == "TRUE" else False
                pattern = re.compile("^({0})$".format(pattern))
                # if name == "clone":
                #     print "DEBUG", name, description, mss_required, entry, value
                field = MetadataField(name, description, qualifier, feature, entry, type_, mss_required, dfast_required, pattern, value)
                fields[name] = field
        for key, value in formData.items():
            field = fields.get(key)
            if field:
                # if key == "clone":
                #     print "DEBUG2", key, value, field.name, field.value, field.pattern.pattern
                field.setValue(value.strip())
                self.fields[key] = field
                # print("DEBUG3", key, value, field.name, field.value, field.pattern.pattern)

    @staticmethod
    def readFromTSV(fileName):

        D = {}
        for line in open(fileName):
            key, value = line.strip("\n").split("\t")
            if key in ["seqTypes", "seqTopologies", "seqNames"]:
                value = value.replace(",", ";")   # for compatibility of DFAST old version
            # field = metadata.fields.get(key)
            if value:
                if key == "source_country":
                    key = "geo_loc_name"  # compatibility for deprecation of "country" in the source feature

                D[key] = value

        return Metadata.initializeFromFormData(D)

    @staticmethod
    def initializeFromFormData(FormData, removePrivate=False):

        # initialize metadata
        metadata = Metadata.initialize()

        # read values from FormData Dictionary
        for key, value in FormData.items():
            field = metadata.fields.get(key)
            if field:
                if field.feature not in ["REFERENCE", "COMMENT"]:
                    field.setValue(value)
            else:
                # skip processing
                continue
                # raise AssertionError("NO VALUE exist in DEFINITION. {}".format(key))
        metadata.addComments(FormData)
        metadata.addReferences(FormData)
        metadata.addExtendedSourceQualifiers(FormData)
        metadata.setDefaults(removePrivate)
        return metadata

    def setDefaults(self, removePrivate=False):
        def removeFields(name):
            field = self.fields.get(name)
            if field:
                del self.fields[name]
            else:
                # print()
                raise IndexError("Field name not found in " + name)

        def setAttribute(name, attr, value):

            field = self.fields.get(name)
            if field is None:
                raise IndexError("Field name not found")
            if attr == "value":
                field.setValue(value)
            else:
                value = True if value == "TRUE" else value
                value = False if value == "FALSE" else value
                if attr == "pattern":
                    value = re.compile(r"^({0})$".format(value))
                # if attr == "mss_required":
                    # print name, attr, value

                setattr(field, attr, value)

        for line in open(self.METADATA_RULES_FILE):
            if line.startswith("#"):
                continue
            cols = line.strip("\n").split("\t")
            name, flag, pattern, target_name, action, value, private = cols[0], cols[1], cols[2], cols[3], cols[4], cols[5], cols[6]
            # print(line),
            private = private == "private"
            field = self.fields.get(name)
            pattern = re.compile(r"^({0})$".format(pattern))
            currentValue = self.getValue(name)
            # print(name, currentValue, pattern.pattern, target_name, value)
            if field and pattern.match(currentValue):

                if private and removePrivate:
                    removeFields(target_name)
                elif action == "DELETE":
                    removeFields(target_name)
                elif action == "SETDEFAULT":
                    self.setDefaultValue(target_name, value)
                elif action == "ADDFIELD":
                    self.addField(target_name, value=None)
                else:
                    setAttribute(target_name, action, value)

            if flag=="reference":
                for i in range(2, self.refNumber + 1):
                    refName = name + ":" + str(i)
                    currentValue = self.getValue(refName)
                    # print(currentValue, refName, target_name)
                    if pattern.match(currentValue):
                        setAttribute(target_name + ":" + str(i), action, value)


            


    # def setDefaults2a(self):

    #     # for compatibility
    #     organism = self.getValue("organism")
    #     if organism:
    #         words = organism.split()
    #         genus = "Genus" if len(words) == 0 else words[0]
    #         species = "sp." if len(words) <= 1 else " ".join(words[1:])
    #         self.setValue("genus", genus)
    #         self.setValue("species", species)
    #     else:
    #         genus = self.getValue("genus", "Genus")
    #         species = self.getValue("species", "sp.")
    #         self.setValue("organism", genus + " " + species)

    #     default_genomes = ["Assembly Method", "Genome Coverage", "Sequencing Technology", "tagset_id", "mol_type", "organism", "strain"]
    #     default_transcriptome = ["division", "keyword", "Assembly Method", "Sequencing Technology", "tagset_id", "mol_type"]
    #     projectType = self.getValue("projectType", "wgs")

    #     if projectType == "gnm":

    #         requiredFields = ["seqTopologies", "seqNames", "seqTypes"] + default_genomes
    #         self.setValue("seqStatus", "complete")
    #         self.setValue("division", "")
    #         self.setValue("keyword", "")
    #         self.setValue("dataType", "")
    #         self.setValue("seqRank", "")

    #         self.setValue("tagset_id", "Genome-Assembly-Data")
    #         self.setValue("mol_type", "genomic DNA")
    #         self.setDefaultValue("locusTagPrefix", "LOCUS")

    #     elif projectType == "wgs":
    #         requiredFields = ["keyword", "type"] + default_genomes
    #         self.setValue("seqStatus", "draft")
    #         self.setValue("division", "")
    #         self.setDefaultValue("keyword", "WGS; STANDARD_DRAFT")
    #         self.setValue("dataType", "WGS")
    #         self.setDefaultValue("seqNames", "Sequence")
    #         self.setValue("seqTopologies", "")
    #         self.setValue("seqTypes", "")
    #         self.setDefaultValue("seqRank", "contig")

    #         self.setValue("tagset_id", "Genome-Assembly-Data")
    #         self.setValue("mol_type", "genomic DNA")
    #         self.setDefaultValue("locusTagPrefix", "LOCUS")
    #         # dynamically set required fields

    #     elif projectType == "tsa":

    #         self.setValue("tagset_id", "Assembly-Data")
    #         self.setValue("mol_type", "transcribed RNA")
    #         self.fields["coverage"].qualifier = "Coverage"
    #         self.setValue("division", "TSA")
    #         self.setValue("keyword", "TSA; Transcriptome Shotgun Assembly")
    #         # print "KEYWORD:", self.getValue("keyword")
    #         self.setValue("dataType", "")
    #         requiredFields = [] + default_transcriptome

    #     else:
    #         # others
    #         requiredFields = []

    #     fields = [field for field in self.fields.values() if field.qualifier in requiredFields]
    #     for field in fields:
    #         field.mss_required = True
    #         # print self.getValue("mol_type")

if __name__ == '__main__':
    metadataFlle = "/Users/ytanizaw/project/labrep_dev/jobs/53eeb0b5-edc5-427d-99c2-b213d6c187a0/result/metadata.txt"
    # metadataFlle = "/Users/ytanizaw/Desktop/meta.txt"
    # metadataFlle = "/Users/ytanizaw/project/labrep_dev/jobs/d4e5bccd-ad54-4355-90ee-0ade4cdd9633/result/metadata.txt"

    # metadata = Metadata.initialize()
    metadata = Metadata.readFromTSV(metadataFlle)

    # metadata.fields["division"].setValue("CON")
    # metadata.fields["bioproject"].value = ["PRJDB00000"]
    # metadata.fields["biosample"].value = ["SAMD90000000", "SAMD90000001"]
    # metadata.fields["sra"].setValue("DRR5000;DRR50002")
    # metadata.fields["sra"].appendValue("DRR99999")
    # metadata.fields["submitter"].appendValue("Tanizawa,Y.")
    # metadata.fields["submitter"].appendValue("Ken,Y.")

    # metadata.toKeyValueFormat("dummy")
    # print(metadata.fields["biosample"].value)
    # print(metadata.fields["biosample"].render())
    # print metadata.validateValues()
    D = {"comment": "OK", "comment:1": "COM1", "comment:2": "COM2; COM2-2", "comment:3": "COMMENT3 is here.; COM3-2; COM3-3"}
    # D = {}
    metadata.addComments(D)
    # print("comment added finished")
    D = {"reference": "REF 1", "pubmed:2": 9999, "reference: 3": "REF3", "author: 3": "Tanizawa, Y.; Yamada, K.", "year: 3": "2016"}

    # metadata.addReferences(D)

    D = {"organelle": "organelle name test A", "source_note": "Note line1; Note line2", "strain": "strain extended", "source_country": "japan"}
    metadata.addExtendedSourceQualifiers(D)

    # outputData = metadata.renderSourceFeature(seqNum=-1, seqRecord="None")
    # for line in outputData:
    #     print(line)

    # outputData = metadata.renderSourceFeature(seqNum=1, seqRecord=None)
    # for line in outputData:
    #     print(line)

    metadata.setValue("division", "CON")
    metadata.setValue("author", "Taro;Lo:")
    # metadata.toTSV("/Users/ytanizaw/Desktop/meta.txt")
    metadata.setDefaults()
    print(metadata.fields["sequencingTechnology"].mss_required)
    print(metadata.fields["sequencingTechnology"].value)
    print(metadata.fields["sequencingTechnology"].setValue("DUMMY"))
    print(metadata.fields["assemblyMethod"].setValue("VELVET"))
    print(metadata.fields["keyword"].setValue("TSA;STANDARD_DRAFT"))
    ret = metadata.validateValues(ignoreNull=False)

    outputData = metadata.renderCommonEntry(includeSource=True)
    # debug
    for line in outputData:
        print(line)

    for key, value in ret.items():
        print(key, value)
    # for key, value in metadata.fields.items():
    # print(str(value))
