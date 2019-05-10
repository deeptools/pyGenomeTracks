from . BedTrack import BedTrack


class TADsTrack(BedTrack):
    SUPPORTED_ENDINGS = ['.domain', '.domains' '.tad', '.tads']
    TRACK_TYPE = 'domains'

    def plot(self, ax, chrom_region, start_region, end_region):
        """
        Plots the boundaries as triangles in the given ax.
        """
        from matplotlib.patches import Polygon
        ymax = 0.001
        valid_regions = 0
        chrom_region = self.check_chrom_str_bytes(self.interval_tree, chrom_region)
        if chrom_region not in self.interval_tree:
            orig = chrom_region
            chrom_region = self.change_chrom_names(chrom_region)
            chrom_region = self.check_chrom_str_bytes(self.interval_tree, chrom_region)
            self.log.info('Chromosome name: {} does not exists. Changing'
                          ' name to {}'.format(orig, chrom_region))
            if chrom_region not in self.interval_tree:
                self.log.error("*Error*\nNeither " + orig + " "
                               "nor " + chrom_region + " exits as a chromosome"
                               " name.\n")
                return
        for region in sorted(self.interval_tree[chrom_region][start_region:end_region]):
            """
                   ______ y2
                  ""
                 "  "
                "    "
               "      "_____ y1
            _____________________
               x1 x2 x3

            """
            x1 = region.begin
            x2 = x1 + float(region.end - region.begin) / 2
            x3 = region.end
            y1 = 0
            y2 = (region.end - region.begin)

            rgb, edgecolor = self.get_rgb_and_edge_color(region.data)

            triangle = Polygon([[x1, y1], [x2, y2], [x3, y1]], closed=True,
                               facecolor=rgb, edgecolor=edgecolor, linewidth=self.properties['line width'])
            ax.add_artist(triangle)
            valid_regions += 1

            if y2 > ymax:
                ymax = y2

        if valid_regions == 0:
            self.log.warning("No regions found for section {}.".format(self.properties['section_name']))

        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            ax.set_ylim(ymax, 0)
        else:
            ax.set_ylim(0, ymax)
