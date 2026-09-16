# -*- coding: utf-8 -*-
from .readTabular import ReadTabular
from .utilities import InputError


class ReadGwas(ReadTabular):
    """
    Reads a GWAS file. The expected fields are:
    chromosome, position, name, and pvalue.

    Example:
    gwas = ReadGwas("file.gwas")
    for record in gwas:
        print(record.chromosome, record.position, record.pvalue)
    """
    # Define alias column names for GWAS
    alias = {
        'chromosome': ['chr', 'chrom'],
        'position': ['bp', 'pos', 'base_pair_location'],
        'variant_id': ['variant', 'id', 'rsid', 'rs_id', 'snp', 'rs', 'rsnum', 'marker', 'markername'],
        'pvalue': ['p', 'pval', 'p-value', 'p_value', 'p.value']
    }
    required_fields = ['chromosome', 'position', 'pvalue']
    required_fields_constr = [str, int, float]
    fields = ['chromosome', 'position', 'variant_id', 'pvalue']

    def __init__(self, file_handle, has_header=False):
        """
        :param file_handle: file handle
        """
        self.has_header = has_header
        super().__init__(file_handle)

        # Adjust length and line_number if header
        if has_header:
            self.line_number += 1
            self.length -= 1

    def adjust_fields(self):
        if self.has_header:
            # revert the alias dictionary
            synonymous = {v: k for k in self.alias for v in self.alias[k]}
            self.header = next(self.file_handle)
            header_fields = self.get_line_data(self.header)
            fields = []
            for hf in header_fields:
                hfl = hf.lower()
                fields.append(synonymous.get(hfl, hfl))
        else:
            fields = self.fields
        # Check all needed fields are present
        # And which position
        positions = []
        for req_f in self.required_fields:
            if not req_f in fields:
                raise InputError(f"The header does not contain any of the following column name: {req_f}, {",".join(alias[req_f])} which is required.")
            else:
                positions.append([i for i, v in enumerate(fields) if v == req_f][0])
        self.used_fields = fields
        self.req_f_pos = positions

    def get_record(self, gwas_line):
        """
        Processes each line from a GWAS file and returns a namedtuple object.

        :param gwas_line: a single line from the GWAS file
        :return: Record object
        """
        line_data = self.get_line_data(gwas_line)

        if len(line_data) < len(self.used_fields):
            if has_header:
                msg = f"The number of fields was anticipated from the following header\n{self.header}."
            else:
                msg = "We expect at least 4 fields, corresponding to: chromosome, position, name, pvalue."
            raise InputError(f"Line {self.line_number} does not have {len(self.used_fields)} fields: {gwas_line}{msg}")
        for i, req_f in enumerate(self.required_fields):
            try:
                line_data[self.req_f_pos[i]] = \
                    self.required_fields_constr[i](line_data[self.req_f_pos[i]])
            except ValueError as e:
                raise InputError(f"Error parsing line {self.line_number}: {gwas_line}\n{e}")

        return self.Record._make(line_data[:len(self.used_fields)])
