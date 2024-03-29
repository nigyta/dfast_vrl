from dataclasses import dataclass, field

# @dataclass
class DFV_WARNING:
    
    def __init__(self, level, name, message, targets=None):
        self.level, self.name, self.message = level, name, message   
        if targets is None:
            self.targets = [] 
        else:
            self.targets = targets

    # level : str = "WARNING"  # INFO, WARNING, CRITICAL
    # name : str = "name"
    # message : str = "message"
    # targets: list = field(default_factory=list)


    def __repr__(self):
        if self.targets:
            return f"[{self.level}] {self.name}:{self.message} {str(self.targets)}]"
        else:
            return f"[{self.level}] {self.name}:{self.message}"

    def add(self, target):
        self.targets.append(target)

    def __str__(self):
        if self.targets:
            return f"[{self.level}] {self.name}:{self.message} {str(self.targets)}]"
        else:
            return f"[{self.level}] {self.name}:{self.message}"

    def to_tuple(self):
        return (self.level, self.name, self.message, self.targets)


INFO_QUERY_MODIFIED = DFV_WARNING(level="INFO", 
                                     name="INFO_QUERY_MODIFIED", 
                                     message="The query genome has been modified during the preprocessing step.")

INFO_SCAFFOLDING_ENABLED = DFV_WARNING(level="INFO", 
                                     name="INFO_SCAFFOLDING_ENABLED", 
                                     message="Scaffolding is enabled. Contigs has been concatenated with runs of Ns of estimated length.")

TRIM_SLASH = DFV_WARNING(level="INFO", 
                                     name="TRIM_SLASH", 
                                     message="Trimmed slash '//' in the query during the processing step.")

TRIM_TERMINAL_N = DFV_WARNING(level="WARNING", 
                                     name="TRIM_TERMINAL_N", 
                                     message="Trimmed terminal Ns in the query during the processing step.")

SEQUENCE_RENAMED = DFV_WARNING(level="WARNING", 
                                     name="SEQUENCE_RENAMED", 
                                     message="Sequence name is too long. Truncated to 50 letters.")


INFO_DRAFT_GENOME = DFV_WARNING(level="INFO", 
                                     name="INFO_DRAFT_GENOME", 
                                     message="The query genome has been classified as a draft genome.")

INFO_NEARLY_COMPLETE_GENOME = DFV_WARNING(level="INFO", 
                                     name="INFO_NEARLY_COMPLETE_GENOME", 
                                     message="The query genome has been classified as a nearly complete genome.")


VADR_ANNOTATION_FAILED = DFV_WARNING(level="CRITICAL", 
                                     name="VADR_ANNOTATION_FAILED", 
                                     message="No biological feature annotated by VADR. Please check the VADR log file carefully.")

INCOMPLETE_CDS_WARNING = DFV_WARNING(level="INFO", 
                                     name="INCOMPLETE_CDS", 
                                     message="CDS is partial or translation contains ambiguous amino acids. Annotated as misc_feature.")

INCOMPLETE_GENOME_WARNING = DFV_WARNING(level="WARNING", 
                                        name="INCOMPLETE_GENOME", 
                                        message="Incomplete query genome.")


MORE_THAN_ONE_MODEL_USED =  DFV_WARNING(level="CRITICAL", 
                                        name="MORE_THAN_ONE_MODEL_USED", 
                                        message="More than one reference model used in VADR. Please check the results and the intermediate files carefully. This may happen in a draft genome.")

MISSING_FEATURES = DFV_WARNING(level="WARNING", 
                               name="MISSING_FEATURES", 
                               message="Expected biological features are missing in the genome.")

PARTIAL_FEATURES = DFV_WARNING(level="WARNING", 
                                        name="PARTIAL_FEATURES", 
                                        message="Annotated features are partial.")

DUPLICATED_FEATURES = DFV_WARNING(level="WARNING", 
                                        name="DUPLICATED_FEATURES", 
                                        message="Biological features are annotated more than once.")

FRAGMENTED_FEATURES = DFV_WARNING(level="WARNING", 
                                        name="FRAGMENTED_FEATURES", 
                                        message="Biological features are split into two or more fragments.")

INCOMPLETE_GENOME_WARNING = DFV_WARNING(level="WARNING", 
                                        name="INCOMPLETE_GENOME", 
                                        message="Incomplete query genome.")

MISC_FEATURES = DFV_WARNING(level="WARNING", 
                                        name="MISC_FEATURES", 
                                        message="Biological features are annotatated as misc_feature.")

VADR_ANNOTATION_FAILED_WITH_CRITICAL_WARNINGS = DFV_WARNING(level="CRITICAL", 
                                     name="VADR_ANNOTATION_FAILED_WITH_CRITICAL_WARNINGS", 
                                     message="VADR failed with critical warnings. Please check VADR result files carefully.",
                                     targets=["-"])


def create_VADR_warning(fail, alert_desc, alert_detail):
    if fail == "yes":
        level = "CRITICAL"
    else:
        level = "INFO"
    return DFV_WARNING(level=level, name="VADR:" + alert_desc, message=alert_detail)

if __name__ == '__main__':

    print(str(DFV_WARNING))





    print(str(DFV_WARNING))
    print(repr(DFV_WARNING))
    print(BROKEN_CDS_WARNING)
    BROKEN_CDS_WARNING.add("test1")
    BROKEN_CDS_WARNING.add("test2")
    print(BROKEN_CDS_WARNING)